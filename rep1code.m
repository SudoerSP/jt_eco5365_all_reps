% TABLE 2 REPLICATION

clear; clc;

data = readtable('pld2023_dataset_2017.csv','VariableNamingRule','preserve');

% Filter the data for the year 2017
% Note: Since the file we are using already has non-2017 data removed, this is
% technically not necessary 
data = data(data.year == 2017, :);

% Omit countries not included in the original table 2 in the paper
validCodes = {'AUS','BEL','CZE','DNK','ESP','FIN','FRA','DEU','GRC','HUN','IRL','ITA','JPN','KOR','NLD','POL','PRT','SVK','SWE','GBR','USA'};
isValid = ismember(data.countrycode, validCodes);
data = data(isValid,:);

% Extract the relevant columns: countrycode, sector, and PPP_va
requiredVars = {'countrycode', 'sector', 'PPP_va'};
if ~all(ismember(requiredVars, data.Properties.VariableNames))
    error('One or more required columns are missing from the dataset.');
end
data = data(:, requiredVars);

data.countrycode = categorical(data.countrycode);
data.sector = categorical(data.sector);
pivotTable = unstack(data, 'PPP_va', 'sector');

% Extract the list of countries and industry names.
countries = pivotTable.countrycode;
allVars = pivotTable.Properties.VariableNames;
industryNames = allVars(2:end);
prodMatrix = table2array(pivotTable(:, 2:end));

% Normalize industries
usIndex = find(countries == 'USA');
if isempty(usIndex)
    error('USA not found in the countrycode column.');
end
% Retrieve the U.S. productivity values for all industries.
usProd = prodMatrix(usIndex, :);
if any(usProd == 0)
    error('Zero productivity value encountered for USA in one or more industries.');
end
Y1 = prodMatrix ./ usProd;
agrIndex = find(strcmpi(industryNames, 'agr'));
if isempty(agrIndex)
    error('The agr industry (representing Food) was not found in the sector column.');
end
agrValues = Y1(:, agrIndex);
if any(agrValues == 0)
    error('Zero productivity value in the agr industry for one or more countries.');
end
Y2 = Y1 ./ agrValues;

% Convert the final normalized matrix back to a table.
normalizedTable = array2table(Y2, 'VariableNames', industryNames);

normalizedTable.countrycode = cellstr(countries);
normalizedTable = movevars(normalizedTable, 'countrycode', 'Before', 1);

% Display the final table.
disp('Replicated Table 2');
disp(normalizedTable);

% EQUATION 18 REPLICATION

gsTable = readtable('pld2023_dataset_2017.csv','VariableNamingRule','preserve');
oecdData = readtable('OECD_2017.csv','VariableNamingRule','preserve');

% Map OECD's 4-digit industry codes to your (agr, min, man, pu) sectors
industryMap = containers.Map( ...
{'D01T03','D05T08','D10T32','D35'}, ...
    {'agr','min','man','pu'} ...
);

% These are the only overlapping sectors
allowedSectors = {'agr','min','man','pu'};

% Keep only rows in gsTable that match an OECD country and an allowed sector
oecdCountriesList = unique(oecdData.COU);
isValidRow = ismember(gsTable.countrycode, oecdCountriesList) ...
          & ismember(gsTable.sector, allowedSectors);
gsTable = gsTable(isValidRow,:);

% Rename columns in oecdData
oecdData.countrycode = oecdData.COU;
oecdData.sector = oecdData.IND;
oecdData.COU = [];
oecdData.IND = [];

% Map OECD's sector codes to our 4-sector classification
rowsToDelete = false(height(oecdData),1);
for i = 1:height(oecdData)
    currSector = oecdData.sector{i};
    if isKey(industryMap, currSector)
        oecdData.sector{i} = industryMap(currSector);
    else
        rowsToDelete(i) = true;
    end
end
oecdData = oecdData(~rowsToDelete,:);

% Expanded table
uniqueCountries = unique(gsTable.countrycode);
nCountries = numel(uniqueCountries);
wideTable = table(uniqueCountries, 'VariableNames', {'countrycode'});

for s = 1:length(allowedSectors)
    sec = allowedSectors{s};
    
    obsProd = nan(nCountries,1);
    gOutput = nan(nCountries,1);
    iv_var = nan(nCountries,1); 

    for j = 1:nCountries
        ctry = uniqueCountries{j};
        idx = strcmp(gsTable.countrycode, ctry) & strcmp(gsTable.sector, sec);
        if sum(idx)==1
            rowData = gsTable(idx,:);
            if rowData.PPP_va > 0
                obsProd(j) = 1 / rowData.PPP_va;
            end
            
            if rowData.xr>0
                % Using the PPP identities to solve for nominal gross output, Y, in LCU:
                numerator = ( (rowData.PPP_y / rowData.PPP_va) - (rowData.PPP_y / rowData.PPP_z ) ) * rowData.VA;
                denominator = 1 - ( rowData.PPP_y / rowData.PPP_z );
                if denominator ~= 0
                    Y_nominal_LCU = numerator / denominator; 
                    % Then convert nominal LCU to thousands of USD:
                    gOutput(j) = (Y_nominal_LCU / rowData.xr) * 1000;
                else
                    gOutput(j) = NaN; 
                end

            end
            
            if rowData.VA>0
                iv_var(j) = log(rowData.VA);
            end
        end
    end
    wideTable.(['observed_prod_' sec]) = obsProd;
    wideTable.(['gross_output_'  sec]) = gOutput;
    wideTable.(['iv_lnVA_' sec]) = iv_var;
end

% Compute domestic share = 1 - import_penetration
% and import_penetration = totalImports / (gross_output + imports - exports)

for s = 1:length(allowedSectors)
    sec = allowedSectors{s};
    
    domShare = nan(nCountries,1);
    
    for j = 1:nCountries
        ctry = uniqueCountries{j};
        
        % Sum of total imports
        idxImp = strcmp(oecdData.countrycode, ctry) ...
               & strcmp(oecdData.sector, sec) ...
               & strcmp(oecdData.Flow, 'Imports');
        totalImports = sum(oecdData.OBS_VALUE(idxImp));
        
        % Sum of total exports
        idxExp = strcmp(oecdData.countrycode, ctry) ...
               & strcmp(oecdData.sector, sec) ...
               & strcmp(oecdData.Flow, 'Exports');
        totalExports = sum(oecdData.OBS_VALUE(idxExp));
        
        % Retrieve the country's gross_output_{sec}
        colGO = ['gross_output_' sec];
        GO_value = wideTable.(colGO)(j);
        
        denom = GO_value + totalImports - totalExports;
        if denom>0
            IPR = totalImports / denom;  % import-penetration ratio
            domShare(j) = 1 - IPR;
        end
    end
    
    wideTable.(['domestic_share_' sec]) = domShare;
end

% Remove rows with any NaN
allVars = wideTable.Properties.VariableNames;
ignoreVars = {'countrycode'};
checkVars = setdiff(allVars, ignoreVars);
rowsBad = false(height(wideTable),1);
for v = 1:length(checkVars)
    colVals = wideTable.(checkVars{v});
    rowsBad = rowsBad | isnan(colVals);
end
wideTable = wideTable(~rowsBad,:);

% Combined table
exportsData = oecdData(strcmp(oecdData.Flow, 'Exports'), :);
exportsData = innerjoin(exportsData, wideTable, ...
    'LeftKeys','countrycode','RightKeys','countrycode');
nObs = height(exportsData);
obsProdVals = nan(nObs,1);
domShareVals = nan(nObs,1);
ivVals = nan(nObs,1);

for i = 1:nObs
    if iscell(exportsData.sector)
        currSec = exportsData.sector{i};
    else
        currSec = char(exportsData.sector(i));
    end
    colObs = ['observed_prod_' currSec];
    colDom = ['domestic_share_' currSec];
    colIV = ['iv_lnVA_' currSec];
    
    obsProdVals(i) = exportsData.(colObs)(i);
    domShareVals(i) = exportsData.(colDom)(i);
    ivVals(i) = exportsData.(colIV)(i);
end

exportsData.observed_prod = obsProdVals;
exportsData.domestic_share = domShareVals;
exportsData.iv_lnVA = ivVals;

exportsData.x_corrected = log(exportsData.OBS_VALUE) - log(exportsData.domestic_share);
exportsData.x_exports = log(exportsData.OBS_VALUE);

exportsData.z = log(exportsData.observed_prod);

% Remove invalid
isValid = true(height(exportsData),1);
isValid(exportsData.OBS_VALUE <= 0) = false;
isValid(exportsData.domestic_share <= 0) = false;
isValid(exportsData.observed_prod <= 0) = false;
badLog = @(v)(isnan(v) | isinf(v));
tmp = badLog(exportsData.x_corrected) | badLog(exportsData.x_exports) | badLog(exportsData.z);
isValid(tmp) = false;

exportsData = exportsData(isValid,:);
exportsData.FE1 = categorical(strcat(exportsData.countrycode, '_', exportsData.PAR));
exportsData.FE2 = categorical(strcat(exportsData.PAR, '_', exportsData.sector));

% OLS
warning('off','stats:LinearModel:RankDefDesignMat');

mdl_ols_corr = fitlm(exportsData, 'x_corrected ~ z + FE1 + FE2 -1');
mdl_ols_exp  = fitlm(exportsData, 'x_exports   ~ z + FE1 + FE2 -1');

warning('on','stats:LinearModel:RankDefDesignMat');

coef_ols_corr = mdl_ols_corr.Coefficients;
theta_ols_corr = coef_ols_corr{'z','Estimate'};
se_ols_corr    = coef_ols_corr{'z','SE'};

coef_ols_exp = mdl_ols_exp.Coefficients;
theta_ols_exp = coef_ols_exp{'z','Estimate'};
se_ols_exp    = coef_ols_exp{'z','SE'};

% IV approach (two-stage least squares)
% We'll use "iv_lnVA" as instrument due to absence of R&D in newer data
% file
% Stage1: z = alpha + gamma * iv_lnVA + FE1 + FE2 + e
% Stage2: x_corrected (or x_exports) = beta + theta * z_hat + FE1 + FE2 + u

warning('off','stats:LinearModel:RankDefDesignMat');
mdl_1st = fitlm(exportsData, 'z ~ iv_lnVA + FE1 + FE2 -1');
warning('on','stats:LinearModel:RankDefDesignMat');

z_hat = mdl_1st.Fitted;

% 2nd stage
exportsData.zhat = z_hat;
warning('off','stats:LinearModel:RankDefDesignMat');
mdl_iv_corr = fitlm(exportsData, 'x_corrected ~ zhat + FE1 + FE2 -1');
mdl_iv_exp = fitlm(exportsData, 'x_exports   ~ zhat + FE1 + FE2 -1');
warning('on','stats:LinearModel:RankDefDesignMat');

coef_iv_corr = mdl_iv_corr.Coefficients;
theta_iv_corr = coef_iv_corr{'zhat','Estimate'};
se_iv_corr = coef_iv_corr{'zhat','SE'};

coef_iv_exp = mdl_iv_exp.Coefficients;
theta_iv_exp = coef_iv_exp{'zhat','Estimate'};
se_iv_exp = coef_iv_exp{'zhat','SE'};

% Display table similar to table 3

ResultTable = table( ...
    [theta_ols_corr; theta_ols_exp; theta_iv_corr; theta_iv_exp], ...
    [se_ols_corr;    se_ols_exp;   se_iv_corr;    se_iv_exp], ...
    {'OLS'; 'OLS'; 'IV'; 'IV'}, ...
    {'Exp×Imp,Ind×Imp FE'; 'Exp×Imp,Ind×Imp FE'; 'Exp×Imp,Ind×Imp FE'; 'Exp×Imp,Ind×Imp FE'}, ...
    'VariableNames', {'Coefficient','StdErr','Method','FixedEffects'} ...
);

dataCell = table2cell(ResultTable);
transposedData = dataCell';
varNamesOriginal = ResultTable.Properties.VariableNames';
newCell = [varNamesOriginal, transposedData];
newHeaders = [{'Variable'}, {'log(corrected exports)'}, {'log(exports)'}, {'log(corrected exports) '}, {'log(exports) '}];
T_transposed = cell2table(newCell, 'VariableNames', newHeaders);

disp('Replication of equation 18');
disp(T_transposed);

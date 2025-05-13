
clear; clc;

WDI = readtable('WDI.csv','VariableNamingRule','preserve');

varsWithYears = WDI.Properties.VariableNames(contains(WDI.Properties.VariableNames,'[YR'));
WDI_years_list_long = unique(cellfun(@(c) sscanf(c,'%d'), varsWithYears));

iso3_country_list = unique(WDI.("Country Code")(~cellfun(@isempty,WDI.("Country Code"))));

nT = numel(WDI_years_list_long);
nC = numel(iso3_country_list);
N = nT * nC;

endowmentFacTable = table( ...
      repmat({''},N,1), ... % iso3
      nan(N,1), ... % year
      nan(N,1), ... % labor
      nan(N,1), ... % land
      nan(N,1), ... % capital_formation
      nan(N,1), ... % gdp
      nan(N,1), ... % capital
      'VariableNames',{'iso3','year','labor','land','capital_formation','gdp','capital'} );

row = 1;
for y = 1:nT
    yval = WDI_years_list_long(y);
    ycol = sprintf('%d [YR%d]',yval,yval);
    for c = 1:nC
        iso = iso3_country_list{c};
        endowmentFacTable.iso3{row} = iso;
        endowmentFacTable.year(row) = yval;

        maskC = strcmp(WDI.("Country Code"),iso);

        m = maskC & strcmp(WDI.("Series Name"),'Labor force, total');
        if any(m); endowmentFacTable.labor(row) = WDI{m,ycol}; end

        m = maskC & strcmp(WDI.("Series Name"),'Agricultural land (sq. km)');
        if any(m); endowmentFacTable.land(row) = WDI{m,ycol}; end

        m = maskC & strcmp(WDI.("Series Name"),'Gross fixed capital formation (constant 2015 US$)');
        if any(m); endowmentFacTable.capital_formation(row) = WDI{m,ycol}; end

        m = maskC & strcmp(WDI.("Series Name"),'GDP (constant 2015 US$)');
        if any(m); endowmentFacTable.gdp(row) = WDI{m,ycol}; end

        row = row + 1;
    end
end

% Capital stock by perpetual‑inventory (5% deprecation)
for c = 1:nC
    iso = iso3_country_list{c};
    rows = find(strcmp(endowmentFacTable.iso3,iso));
    [~,o]= sort(endowmentFacTable.year(rows)); rows = rows(o);

    x = endowmentFacTable.capital_formation(rows);

    g = 0;
    if numel(x) >= 6 && all(~isnan(x(1:6))) && all(x(1:5)>0)
        gVec = diff(x(1:6)) ./ x(1:5);
        g = mean(gVec,'omitnan');
    end

    x0 = x(1);
    if g > 0
        K0 = (1+g)/g * x0;
    else
        % fallback: constant flow
        K0 = x0;
    end

    K = nan(size(x));
    K(1) = K0;
    for t = 2:numel(x)
        K(t) = 0.95*K(t-1) + x(t);
    end
    endowmentFacTable.capital(rows) = K;
end

endowmentFacTable(cellfun(@isempty,endowmentFacTable.iso3),:) = [];
warning('off', 'MATLAB:nearlySingularMatrix');

oecdFiles = {'OECD_AGRI.csv','OECD_ELEN.csv','OECD_MINING.csv', ...
             'OECD_PETROL.csv','OECD_TEXT.csv','OECD_TRANSPORT.csv', ...
             'OECD_WOOD.csv'};
usFiles = {'US_AGRI.csv','US_ELEN.csv','US_MINING.csv','US_PETROL.csv', ...
             'US_TEXT.csv','US_TRANSPORT.csv','US_WOOD.csv'};
indLbl = {'AGRI','ELEN','MINING','PETROL','TEXT','TRANSPORT','WOOD'};

if ismember('gpd',endowmentFacTable.Properties.VariableNames)
    endowmentFacTable.Properties.VariableNames{'gpd'} = 'gdp';
end
years_Endow = unique(endowmentFacTable.year);

Tg = readtable('CEPII_Gravity_Filtered.csv','VariableNamingRule','preserve');
years_Grav = unique(Tg.year);

years_OECD = years_Grav;
for f = 1:numel(oecdFiles)
    To = readtable(oecdFiles{f},'VariableNamingRule','preserve');
    years_OECD = intersect(years_OECD,unique(To.TIME_PERIOD));
end

years_US = years_Grav;
for f = 1:numel(usFiles)
    Tu = readtable(usFiles{f},'VariableNamingRule','preserve');
    years_US = intersect(years_US,unique(Tu.Time));
end

commonYears = sort(intersect(intersect(years_Endow,years_Grav), ...
                             intersect(years_OECD ,years_US  )));

T_VA = [];
for k = 1:numel(oecdFiles)
    To = readtable(oecdFiles{k},'VariableNamingRule','preserve');
    To = To(ismember(To.TIME_PERIOD,commonYears) & strcmp(To.Measure,'Value added'),:);
    To.iso3 = string(To.REF_AREA); To.year = To.TIME_PERIOD;
    vName = ['ValueAdded_' indLbl{k}]; To.(vName) = To.OBS_VALUE;
    To = groupsummary(To(:,{'iso3','year',vName}),{'iso3','year'},'sum',vName);
    To.Properties.VariableNames(end) = {vName};
    if k==1, T_VA = To; else, T_VA = innerjoin(T_VA,To,'Keys',{'iso3','year'}); end
end

selfIdx = strcmp(Tg.iso3_o,Tg.iso3_d) & ismember(Tg.year,commonYears);
Tland = Tg(selfIdx,{'iso3_o','distcap'});
Tland.land = pi*((3/2)*Tland.distcap).^2;
Tland = groupsummary(Tland,'iso3_o','mean','land');
Tland.Properties.VariableNames = {'iso3','GroupCount','land'}; Tland.GroupCount = [];

T_end = endowmentFacTable(ismember(endowmentFacTable.year,commonYears),:);
T_end = rmmissing(T_end,'DataVariables',{'capital','labor','gdp'});
T_end = removevars(T_end,'land');
T_end = innerjoin(T_end,Tland,'Keys','iso3');

T_end.ratioCapitalLand = T_end.capital ./ T_end.land;
T_end.ratioLaborLand = T_end.labor ./ T_end.land;
T_end.logCapitalLand = log(T_end.ratioCapitalLand);
T_end.logLaborLand = log(T_end.ratioLaborLand);
T_end = rmmissing(T_end);

Tdist = Tg(strcmp(Tg.iso3_o,'USA') & ~strcmp(Tg.iso3_d,'USA') & ...
           ismember(Tg.year,commonYears),{'iso3_d','year','distcap'});
Tdist.Properties.VariableNames = {'iso3','year','dist_to_US'};

map = { ...
    'Australia', 'AUS'; 'Austria', 'AUT'; 'Belgium', 'BEL'; ...
    'Canada', 'CAN'; 'Chile', 'CHL'; 'Colombia', 'COL'; ...
    'Costa Rica', 'CRI'; 'Czech Republic', 'CZE'; 'Denmark', 'DNK'; ...
    'Estonia', 'EST'; 'Finland', 'FIN'; 'France', 'FRA'; ...
    'Germany', 'DEU'; 'Greece', 'GRC'; 'Hungary', 'HUN'; ...
    'Iceland', 'ISL'; 'Ireland', 'IRL'; 'Israel', 'ISR'; ...
    'Italy', 'ITA'; 'Japan', 'JPN'; 'Korea, South', 'KOR'; ...
    'Latvia', 'LVA'; 'Lithuania', 'LTU'; 'Luxembourg', 'LUX'; ...
    'Mexico', 'MEX'; 'Netherlands', 'NLD'; 'New Zealand', 'NZL'; ...
    'Norway', 'NOR'; 'Poland', 'POL'; 'Portugal', 'PRT'; ...
    'Slovakia', 'SVK'; 'Slovenia', 'SVN'; 'Spain', 'ESP'; ...
    'Sweden', 'SWE'; 'Switzerland', 'CHE'; 'Turkey', 'TUR'; ...
    'United Kingdom', 'GBR'; 'United States', 'USA'  };
name2iso = containers.Map(map(:,1),map(:,2));
cleanNum = @(x) str2double(erase(string(x),','));

for f = 1:numel(usFiles)
    Tu = readtable(usFiles{f},'VariableNamingRule','preserve');
    Tu = Tu(ismember(Tu.Time,commonYears) & isKey(name2iso,Tu.Country),:);
    Tu.iso3 = cellfun(@(c) name2iso(c),Tu.Country,'uni',false); Tu.year = Tu.Time;
    Tu.val = cleanNum(Tu.("Customs Value (Gen) ($US)"));

    grp = groupsummary(Tu,{'iso3','year'},'sum','val');
    grp.Properties.VariableNames{'sum_val'} = 'expTot';
    world = groupsummary(Tu,'year','sum','val'); world.Properties.VariableNames{'sum_val'}='worldTot';
    grp = innerjoin(grp,world,'Keys','year');
    grp.logVar = log(grp.expTot ./ grp.worldTot);

    vName = ['LogExportVar_' indLbl{f}];
    grp = grp(:,{'iso3','year','logVar'}); grp.Properties.VariableNames{'logVar'} = vName;

    if f==1, T_var = grp; else, T_var = innerjoin(T_var,grp,'Keys',{'iso3','year'}); end
end

T = innerjoin(T_VA,T_var,'Keys',{'iso3','year'});
T = innerjoin(T,T_end ,'Keys',{'iso3','year'});
T = innerjoin(T,Tdist ,'Keys',{'iso3','year'});

vaNames = strcat('ValueAdded_',indLbl);
T.TotalValueAdded = sum(T{:,vaNames},2,'omitnan');
for k = 1:numel(indLbl)
    T.( ['Share_' indLbl{k}] ) = T.(vaNames{k}) ./ T.TotalValueAdded;
end

ordered = [ {'iso3','year'}, vaNames, {'TotalValueAdded'}, strcat('Share_',indLbl), ...
            strcat('LogExportVar_',indLbl), ...
            {'land','capital','labor','gdp','ratioCapitalLand','ratioLaborLand', ...
             'logCapitalLand','logLaborLand','dist_to_US'} ];

T_translog = T(:,ordered); T_translog = rmmissing(T_translog);

sectorLabelsForStudy = {'AGRI','ELEN','MINING','PETROL','TEXT','TRANSPORT','WOOD'};
countOfSectors = numel(sectorLabelsForStudy);
% 7 share and 1 TFP
countOfEquations = countOfSectors + 1;

sharesMatrix =  T_translog{:, strcat('Share_', sectorLabelsForStudy)};
exportVarMatrix =  T_translog{:, strcat('LogExportVar_',sectorLabelsForStudy)};
distanceInKmToUSA =  T_translog.dist_to_US;
distanceSquared =  distanceInKmToUSA.^2;
logLaborPerLand =  T_translog.logLaborLand;
logCapitalPerLand =  T_translog.logCapitalLand;
logGDPperWorker =  log(T_translog.gdp ./ T_translog.labor);
rowCountTotal =  height(T_translog);

% common instruments that are available in‑memory (constant included)
Z_common = [ ...
        ones(rowCountTotal,1), ...
        distanceInKmToUSA, ...
        distanceSquared, ...
        logLaborPerLand, ...
        logCapitalPerLand ];

% First stage
fittedExportVarMatrix = zeros(rowCountTotal, countOfSectors);

for idxSector = 1:countOfSectors
    y_variety = exportVarMatrix(:,idxSector);
    % OLS β̂_zi
    beta_Z = (Z_common' * Z_common) \ (Z_common' * y_variety);
    fittedExportVarMatrix(:,idxSector) = Z_common * beta_Z;
    % ŷ = Zβ̂_zi
end

bigCell_X = cell(countOfEquations,1);
bigCell_y = cell(countOfEquations,1);
residualMatrix = zeros(rowCountTotal, countOfEquations);

% share equations
for eqIdx = 1:countOfSectors
    % dependent variable
    bigCell_y{eqIdx} = sharesMatrix(:,eqIdx);

    % regressors: constant | seven fitted logs | exogenous instruments
    regressorBlock = [ones(rowCountTotal,1), fittedExportVarMatrix, Z_common(:,2:end)];
    bigCell_X{eqIdx} = regressorBlock;

    % 2‑SLS coefficients (X already orthogonal to the instrumented vars)
    betaTS = (regressorBlock' * regressorBlock) \ ...
        (regressorBlock' * bigCell_y{eqIdx});

    % residuals for Ω̂
    residualMatrix(:,eqIdx)= bigCell_y{eqIdx} - regressorBlock*betaTS;
end

% the TFP equation
bigCell_y{countOfEquations} = logGDPperWorker;
bigCell_X{countOfEquations} = [ones(rowCountTotal,1), fittedExportVarMatrix, Z_common(:,2:end)];

betaTS_TFP = (bigCell_X{end}'*bigCell_X{end}) \ ...
    (bigCell_X{end}'*logGDPperWorker);
residualMatrix(:,end) = logGDPperWorker - bigCell_X{end}*betaTS_TFP;

% covariance of the 2SLS residuals (Ω̂)
omegaHat = (residualMatrix' * residualMatrix) ./ rowCountTotal;

% full 3‑SLS
inverseOmegaHat = inv(omegaHat);

dimensionTotals = cellfun(@(x) size(x,2), bigCell_X);
parameterPositions = [0; cumsum(dimensionTotals(:))];
grandK = parameterPositions(end);

XX_weighted = zeros(grandK, grandK);
Xy_weighted = zeros(grandK, 1);

for i = 1:countOfEquations
    Xi = bigCell_X{i};
    yi = bigCell_y{i};
    idxI = (parameterPositions(i)+1) : parameterPositions(i+1);

    for j = 1:countOfEquations
        Xj = bigCell_X{j};
        idxJ = (parameterPositions(j)+1) : parameterPositions(j+1);
        w_ij = inverseOmegaHat(i,j);

        XX_weighted(idxI,idxJ) = XX_weighted(idxI,idxJ) + w_ij * (Xi' * Xj);
        Xy_weighted(idxI) = Xy_weighted(idxI) + w_ij * (Xi' * bigCell_y{j});
    end
end

betaThreeSLS = XX_weighted \ Xy_weighted;


% extract γ̂ AND ρ̂ and impose constraints

offsetLogs = 1;
gammaMatrixHat = zeros(countOfSectors);
for i = 1:countOfSectors
    eqStart = parameterPositions(i) + 1;
    gammaMatrixHat(i,:)= betaThreeSLS(eqStart+offsetLogs : ...
                                      eqStart+offsetLogs+countOfSectors-1).';
end

% symmetry (γij = γji except on diagonals)
gammaMatrixHat = 0.5*(gammaMatrixHat + gammaMatrixHat.');

gammaMatrixHat(end,:) = -sum(gammaMatrixHat(1:end-1,:),1);

% ρ̂_i
eqStart_TFP = parameterPositions(end-1) + 1;
rhoVectorHat = betaThreeSLS(eqStart_TFP+offsetLogs : ...
                                eqStart_TFP+offsetLogs+countOfSectors-1);

tableOfReplication = array2table([gammaMatrixHat, rhoVectorHat], ...
                         'VariableNames', ...
                         [sectorLabelsForStudy, {'TFP'}], ...
                         'RowNames', sectorLabelsForStudy);

disp('================================ Table 3 ================================');
disp(tableOfReplication);

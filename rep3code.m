
clear;  clc;

opts = {'VariableNamingRule','preserve'};
comtrade = readtable('Comtrade_2023.csv', opts{:});
dist     = readtable('dist_cepii.csv'  , opts{:});
wdi      = readtable('WDI.csv'         , opts{:});

comtrade = comtrade(strcmp(comtrade.flowDesc,'Export'), :);
need = {'reporterISO','partnerISO','primaryValue'};
for k = 1:numel(need)
    comtrade = comtrade(~ismissing(comtrade.(need{k})), :);
end
comtrade = comtrade(comtrade.primaryValue >= 0, :);

dist = dist(~ismissing(dist.iso_o) & ~ismissing(dist.iso_d) ...
            & ~isnan(dist.dist)    & ~isnan(dist.contig), :);

wdi = wdi(strcmp(wdi.('Series Code'),'NY.GDP.MKTP.CD'), :);
wdi = wdi(~isnan(wdi.('2023 [YR2023]')), :);

codesT = unique([comtrade.reporterISO ; comtrade.partnerISO]);
codesD = unique([dist.iso_o ; dist.iso_d]);
codesG = unique(wdi.('Country Code'));
iso3   = intersect(intersect(codesT,codesD), codesG);

fprintf('ISO-3 codes present in *all* datasets (%d countries):\n',numel(iso3));
fprintf('%s\n', strjoin(iso3, ', '));

inSet   = @(x) ismember(x,iso3);
comtrade = comtrade( inSet(comtrade.reporterISO) & inSet(comtrade.partnerISO), :);
dist     = dist(     inSet(dist.iso_o)           & inSet(dist.iso_d)         , :);
wdi      = wdi(      inSet(wdi.('Country Code')), :);

[G,eISO,iISO] = findgroups(comtrade.reporterISO, comtrade.partnerISO);
pairVal       = splitapply(@sum, comtrade.primaryValue, G);
pairs         = table(eISO,iISO,pairVal, ...
                      'VariableNames',{'expISO','impISO','exportsUSD'});

gdpMap   = containers.Map(wdi.('Country Code'), wdi.('2023 [YR2023]'));
distKey  = strcat(dist.iso_o,'_',dist.iso_d);
distMap  = containers.Map(distKey , dist.dist);
bordMap  = containers.Map(distKey , dist.contig);

worldGDP = sum(cell2mat(values(gdpMap)));

remMap = containers.Map('KeyType','char','ValueType','double');
for c = 1:numel(iso3)
    i = iso3{c};
    s = 0;
    for d = 1:numel(iso3)
        j = iso3{d};
        if strcmp(i,j), continue; end
        key = sprintf('%s_%s',i,j);
        if ~isKey(distMap,key), key = sprintf('%s_%s',j,i); end
        if ~isKey(distMap,key), continue; end
        s = s + distMap(key) * (gdpMap(j)/worldGDP);
    end
    remMap(i) = max(s,1);
end

n = height(pairs);
vars = {'ln_exports','ln_gdp_exp','ln_gdp_imp', ...
        'ln_dist','ln_rem_exp','ln_rem_imp','border'};
T = array2table(NaN(n,numel(vars)), 'VariableNames',vars);

for r = 1:n
    ei = pairs.expISO{r};   ii = pairs.impISO{r};
    % GDP
    g_e = gdpMap(ei);   g_i = gdpMap(ii);
    % Distance + border
    key = sprintf('%s_%s',ei,ii);
    if ~isKey(distMap,key), key = sprintf('%s_%s',ii,ei); end
    d    = distMap(key);
    bdum = bordMap(key);
    % Remoteness
    r_e = remMap(ei);   r_i = remMap(ii);
    % Logs (+1 USD keeps zero flows)
    T.ln_exports(r) = log(pairs.exportsUSD(r) + 1);
    T.ln_gdp_exp(r) = log(g_e);
    T.ln_gdp_imp(r) = log(g_i);
    T.ln_dist(r)    = log(d);
    T.ln_rem_exp(r) = log(r_e);
    T.ln_rem_imp(r) = log(r_i);
    T.border(r)     = bdum;
end

spec = 'ln_exports ~ ln_gdp_exp + ln_gdp_imp + ln_dist + ln_rem_exp + ln_rem_imp + border';
mdl  = fitlm(T, spec);

disp('OLS results for 2023 cross-section (classic s.e.’s):');
disp(mdl);

coef = mdl.Coefficients;
bExp  = coef.Estimate(strcmp(coef.Properties.RowNames,'ln_gdp_exp'));
bImp  = coef.Estimate(strcmp(coef.Properties.RowNames,'ln_gdp_imp'));
bDist = coef.Estimate(strcmp(coef.Properties.RowNames,'ln_dist'));

fprintf('\n--- Simple diagnostics (descriptive only) ---\n');
fprintf('Exporter-GDP elasticity (coef on ln GDP_exp, expect near 1): %.3f\n', bExp);
fprintf('Importer-GDP elasticity (coef on ln GDP_imp, expect near 1): %.3f\n', bImp);
fprintf('Distance coefficient (coef on ln dist, expect near -1):       %.3f\n', bDist);

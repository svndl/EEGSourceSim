function [data,oltrreport]=outliertreatment(data,options)
%OUTLIERTREATMENT  Treats outliers, mostly with preprocessing
%   [DATA,REPORT] = OUTLIERTREATMENT(DATA,OPTIONS), where DATA is a regular
%   N-by-3 data set and OPTIONS is a valid options structure for DMAT,
%   returns DATA, which may contain fewer data points than before, and
%   REPORT, a structure with information regarding the outlier treatment
%   process.
%
%   See also EWMAV2, PLOTEWMA.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.3: Edited so relative cut-off values are returned to oltrreport.

%% Initialize some variables
dmatewmaplot=[];
co=0;

%% Depending on options, treat outliers
switch options.OutlierTreatment(1),
    case {'R','r'}
        oltrreport.method='Relative cut-off';
        me = mean(data(:,3));
        sd = std(data(:,3));
        options.OutlierMax = me+options.OutlierMax*sd;
        options.OutlierMin = me+options.OutlierMin*sd;
        use = data(:,3)<options.OutlierMax & data(:,3)>options.OutlierMin;
        data = data(use,:);
    case {'A','a'}
        oltrreport.method='Absolute cut-off';
        use = data(:,3)<options.OutlierMax & data(:,3)>options.OutlierMin;
        data = data(use,:);
    case {'M','m'}
        oltrreport.method='Mixed model';
        use = true(length(data),1);
    case {'E','e'}
        oltrreport.method='EWMA';
        [co dmatewmaplot] = ewmav2(data,options.EWMA.L,...
            options.EWMA.l,options.EWMA.s);
        use = data(:,3)>co;
        data = data(use,:);
    case {'B','b'}
        oltrreport.method='EWMA and Mixed model';
        [co dmatewmaplot]= ewmav2(data,options.EWMA.L,...
            options.EWMA.l,options.EWMA.s);
        use = data(:,3)>co;
        data = data(use,:);
    case {'N','n'}
        oltrreport.method='Outlier treatment: None';
        use = true(length(data),1);
end

%% Record how much was censored
conds = unique(data(:,1));
nsets = length(conds);
used = zeros(1,nsets);
for ctr=1:nsets
    used(ctr) = 100*mean(use(data(:,1)==conds(ctr)));
end

%% Prepare output
oltrreport = struct('cutoff',co,...
    'used',used,...
    'use',use,...
    'EWMA',options.EWMA,...
    'EWMAplot',dmatewmaplot,...
    'OutlierMax',options.OutlierMax,...
    'OutlierMin',options.OutlierMin);
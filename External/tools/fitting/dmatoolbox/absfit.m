function fstr = absfit(param,data,k)
%ABSFIT  Calculates model fit measures
%   FSTRUCT = ABSFIT(PARAM,DATA,DF), where PARAM is a full parameter set
%   (usually the estimated parameter set), DATA is a cell matrix where
%   each cell represents a condition, and each cell is formatted as
%   explained in the help for FITDIFFV13, and DF is the number of degrees
%   of freedom in the model, returns FSTRUCT, a structure with the
%   following fields:
%      .IgnoredConditions         Logical index to any conditions that have
%                                 been ignored during the fitting (due to
%                                 too few data points). Column 1 is for X=0
%                                 responses and Column 2 for X=1.
%      .ObservedFrequencies       The observed response frequencies, split
%                                 by condition, bin, and response.
%                                 Odd-numbered rows are for X=0 responses,
%                                 even-numbered for X=1. Rows 1 and 2 are
%                                 for condition 1, 3 and 4 for condition 2
%                                 etc. Columns are for bins.
%      .ObservedProportions       The observed response proportions.
%      .PredictedFrequencies      The predicted response frequencies.
%      .PredictedProportions      The predicted response frequencies.
%      .LogLikelihoodModel        -2*log-likelihood of the data, given the
%                                 model.
%      .LogLikelihoodData         -2*log-likelihood of the data, given the
%                                 saturated model.
%      .LogLikelihoodRatio        The log-likelihood ratio statistic G2.
%      .ChiSquare                 The Chi-Square statistic X2.
%      .AIC                       The Akaike Information Criterion.
%      .AICc                      The Small-Sample Akaike Information
%                                 Criterion.
%      .BIC                       The Bayesian Information Criterion.
%
%   See also MODELFITTABLE and MULTIESTV4.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.4: Fixed formulae.

nsets = size(param,1);
nbins = length(data{1}{1})+1;

c = [data{:}];
ns = [c{7,:}];
N = sum(ns(:));
nprop = repmat(reshape([ns(:) ns(:)]',2*nsets,1)./N,1,nbins);
u = false(2*nsets,1);
predquant = zeros(2*nsets,5); perc = .1:.2:.9;
allPp = zeros(nsets*2,nbins);
allOn = zeros(nsets*2,nbins);
usedN = N;
for ctr=1:nsets
    cc = 2*ctr+[-1 0];
    [Pcond u(cc) deleted] = condition_n(param(ctr,:),data{ctr},nbins);
    usedN = usedN-deleted;
    allPp(cc,:) = Pcond;
    if isempty(data{ctr}{3}), data{ctr}{3}=zeros(1,nbins); end
    if isempty(data{ctr}{4}), data{ctr}{4}=zeros(1,nbins); end
    allOn(cc,:) = [data{ctr}{4}(:)';data{ctr}{3}(:)'];
    predquant(cc,:) = [invcdfdif(perc,0,param(ctr,:),1e-8);
        invcdfdif(perc,1,param(ctr,:),1e-8)];
end

allPn = allPp.*usedN.*nprop;
allPn(~u,:) = NaN;
allOn(~u,:) = NaN;
allOp = allOn./(usedN.*nprop);
allOp(~u,:) = NaN;
allPp(~u,:) = NaN;

allOp = allOn./(usedN.*nprop);
Nvalue = 1e-5;
allPp(allPp<Nvalue)=Nvalue;
allOp(allOp<Nvalue)=Nvalue;
% allPn = allPp.*(usedN.*nprop);
% allOn = allOp.*(usedN.*nprop);

xn = reshape(allOn(u,:),[],1);
mn = reshape(allPn(u,:),[],1);
xp = reshape(allOp(u,:),[],1);
mp = reshape(allPp(u,:),[],1);

Lm = sum(xn.*log(mp));
Lx = sum(xn.*log(xp));
LR = -2*sum(xn.*log(mp./xp));

fstr.IgnoredConditions = reshape(~u,2,[])';
fstr.ObservedFrequencies = allOn;
fstr.ObservedProportions = allOp;
fstr.PredictedFrequencies = allPn;
fstr.PredictedProportions = allPp;
fstr.PredictedQuantiles = predquant;
fstr.LogLikelihoodModel = Lm;
fstr.LogLikelihoodData = Lx;
fstr.LogLikelihoodRatio = LR;
fstr.Cells = numel(xn);
fstr.ChiSquare = sum((xn-mn).^2./mn);
fstr.AIC = -2*Lm + 2*k;
fstr.AICc = -2*Lm + 2*k*N/(N-k-1);
fstr.BIC = -2*Lm + k*log(N);

%% Find expected proportions for each condition
function [Pcond use deleted] = condition_n(param,data,nb)
use = true(2,1);
deleted = 0;
nneeded = 11;
if data{6}<=nneeded 
    Pn = zeros(1,nb); % if number of zeros <=11, ignore
    use(1) = false;
    deleted = deleted + data{6};
else
    Pn = response_n(param,data{2},0,data{8});
end
if data{5}<=nneeded 
    Py = zeros(1,nb); % if number of ones <=11, ignore
    use(2) = false;
    deleted = deleted + data{5};
else
    Py = response_n(param,data{1},1,data{8});
end
Pcond = [Pn;Py];


%% Find expected proportions for each response
function P = response_n(param,Q,reply,extrt)
[P px] = cdfdif(Q,reply,param); % param([8:end]) are ignored
px = 2*reply*px-reply-px+1;
if length(param)>7 && param(8)<1 % If requested, apply mixed model
    olpi = param(8);
    gamm = param(9); % You can't enter 8 parameters
    G = olpi.*[P px]+(1-olpi).*[(Q-extrt(1))./(extrt(2)-extrt(1)) 1].*...
        (gamm/2+(1-gamm).*px);
    px = G(end);
    P = G(1:end-1);
end
P = diff([0 P px]);
P(P<1e-5)=1e-5;

%{
%% Obsolete file:
function fstr = absfit(param,data,k)
%ABSFIT  Calculates model fit measures
%   FSTRUCT = ABSFIT(PARAM,DATA,DF), where PARAM is a full parameter set
%   (usually the estimated parameter set), DATA is a cell matrix where
%   each cell represents a condition, and each cell is formatted as
%   explained in the help for FITDIFFV13, and DF is the number of degrees
%   of freedom in the model, returns FSTRUCT, a structure with the
%   following fields:
%      .IgnoredConditions         Logical index to any conditions that have
%                                 been ignored during the fitting (due to
%                                 too few data points). Column 1 is for X=0
%                                 responses and Column 2 for X=1.
%      .ObservedFrequencies       The observed response frequencies, split
%                                 by condition, bin, and response.
%                                 Odd-numbered rows are for X=0 responses,
%                                 even-numbered for X=1. Rows 1 and 2 are
%                                 for condition 1, 3 and 4 for condition 2
%                                 etc. Columns are for bins.
%      .ObservedProportions       The observed response proportions.
%      .PredictedFrequencies      The predicted response frequencies.
%      .PredictedProportions      The predicted response frequencies.
%      .LogLikelihoodModel        -2*log-likelihood of the data, given the
%                                 model.
%      .LogLikelihoodData         -2*log-likelihood of the data, given the
%                                 saturated model.
%      .LogLikelihoodRatio        The log-likelihood ratio statistic G2.
%      .ChiSquare                 The Chi-Square statistic X2.
%      .AIC                       The Akaike Information Criterion.
%      .AICc                      The Small-Sample Akaike Information
%                                 Criterion.
%      .BIC                       The Bayesian Information Criterion.
%
%   See also MODELFITTABLE and MULTIESTV4.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

nsets = size(param,1);
nbins = max(max(cellfun('length',[data{:}])));

c = [data{:}];
ns = [c{7,:}];
N = sum(ns(:));
nprop = repmat(reshape([ns(:) ns(:)]',2*nsets,1)./N,1,nbins);
u = false(2*nsets,1);
allPp = zeros(nsets*2,nbins);
allOn = zeros(nsets*2,nbins);

for ctr=1:nsets
    cc = 2*ctr+[-1 0];
    [Pcond u(cc)] = condition(param(ctr,:),data{ctr},nbins);
    allPp(cc,:) = Pcond;
    if isempty(data{ctr}{3}), data{ctr}{3}=zeros(1,nbins); end
    if isempty(data{ctr}{4}), data{ctr}{4}=zeros(1,nbins); end
    allOn(cc,:) = [data{ctr}{4};data{ctr}{3}];
end
allPn = allPp.*N.*nprop;
allOp = allOn./sum(allOn(:));
allPn(~u,:) = NaN;
allOp(~u,:) = NaN;
allOn(~u,:) = NaN;
allPp(~u,:) = NaN;

value = .1;
Nvalue = value/N;
allPp(allPp<Nvalue)=Nvalue;
allOp(allOp<Nvalue)=Nvalue;
allPn(allPn<value)=value;
allOn(allOn<value)=value;

xn = reshape(allOn(u,:),[],1);
mn = reshape(allPn(u,:),[],1);

Lm = sum(xn.*log(mn));
Lx = sum(xn.*log(xn));
LR = -2*sum(xn.*log(mn./xn));

fstr.IgnoredConditions = reshape(~u,2,[])';
fstr.ObservedFrequencies = allOn;
fstr.ObservedProportions = allOp;
fstr.PredictedFrequencies = allPn;
fstr.PredictedProportions = allPp;
fstr.LogLikelihoodModel = Lm;
fstr.LogLikelihoodData = Lx;
fstr.LogLikelihoodRatio = LR;
fstr.Cells = numel(xn);
fstr.ChiSquare = sum((xn-mn).^2./mn);
fstr.AIC = Lm+2*k;
fstr.AICc = Lm+2*k*(k+1)/(N-k-1);
fstr.BIC = Lm+k*log(N);


%% Find expected proportions for each condition
function [Pcond use] = condition(param,data,nb)
use = true(2,1);
if data{6}<12
    Pn = zeros(1,nb); % if number of zeros <=11, ignore
    use(1) = false;
else
    Pn = response(param,data{2},0,data{8});
end
if data{5}<12
    Py = zeros(1,nb); % if number of ones <=11, ignore
    use(2) = false;
else
    Py = response(param,data{1},1,data{8});
end
Pcond = [Pn;Py];


%% Find expected proportions for each response
function P = response(param,Q,reply,extrt)
[P px] = cdfdif(Q,reply,param); % param([8:end]) are ignored
px = 2*reply*px-reply-px+1;
if length(param)>7 && param(8)<1 % If requested, apply mixed model
    olpi = param(8);
    gamm = param(9); % You can't enter 8 parameters
    G = olpi.*[P px]+(1-olpi).*[(Q-extrt(1))./(extrt(2)-extrt(1)) 1].*...
        (gamm/2+(1-gamm).*px);
    px = G(end);
    P = G(1:end-1);
end
P = diff([0 P px]);

%}
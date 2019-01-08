function L = multifitv4(small,data,controls,method)
%MULTIFITV4  Compute deviance for all conditions
%   L = MULTIFITV4(SMALL,DATA,CONTROLS,METHOD), computes a deviance for the
%   diffusion model. SMALL is a vector of free parameters, as given by the
%   optimization routine. DATA is a cell matrix where each cell represents
%   a condition, and each cell is formatted as explained in the help for
%   FITDIFFV13. CONTROLS is a structure as returned from function SMALLER.
%   METHOD is a boolean indicating which fit procedure to use. If METHOD is
%   false, a Chi-Square loss is returned. If it is true, a deviance measure
%   based on a multinomial likelihood function (quantile maximum
%   likelihood / quantile probability products) is calculated (-2*log(L)).
%
%   See also FITDIFFV13, SMALLER, BIGGER.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%% Recover parameters
controls.small = small;
controls = bigger(controls);
nsets = size(controls.large,1);

%% Shouldn't be any NaNs or Infs in the parameter set
if any(isinf(controls.large(:))) || any(isnan(controls.large(:)))
    warning('DMAT:multifitv4:NaNInfParameter',...
        'NaN or Inf found in parameter set.')
    L=1e10;
    return
end

%% Is the parameter set possible?
if ~isgood(controls.large)
    L=1e10;
    return
end

%% Fit the model, one condition at a time;
L = .0;
for ctr=1:nsets
    xx = fitdiffv13(controls.large(ctr,:),data{ctr},method);
    if xx<0,xx=1e10; end
    L = L+xx;
end

%% Loss shouldn't be NaN or Inf either
if isnan(L)||isinf(L)||L<0
    warning('DMAT:multifitv4:NaNInfLoss','Loss was %f.',L)
    L=1e10;
end
%fprintf('Objective: %16.6f\n',L)

drawnow % for GUI
function D = fitdiffv13(param,data,method)
%FITDIFFV13  Calculate deviance for one condition
%   D = FITDIFFV13(PARAM,DATA,METHOD), returns a scalar deviance measure.
%
%   PARAM is a possible parameter vector for the present condition. If
%   PARAM has 7 elements, the Ratcliff diffusion model is used, if it has
%   9 elements, a mixed model is used.
%
%   DATA is a cell matrix with statistics for the present data set. DATA
%   has 8 fields: (1) bin edges for "1" responses, (2) bin edges for "0"
%   responses, (3) observed frequencies for "1" responses, (4) observed
%   frequencies for "0" responses, (5) total number of "1" responses, (6)
%   total number of "0" responses, (7) total number of observations, (8) a
%   vector with the lowest and highest response times in the present
%   condition.
%
%   METHOD is a boolean indicating which fit procedure to use. If METHOD
%   is false, a Chi-Square loss is returned. If it is true, a deviance
%   measure based on a multinomial likelihood function (quantile maximum
%   likelihood) is calculated (-2*log(L)).
%
%   See also MULTIFITV4.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Made number of needed responses variable, for future use.

%% Structure of the DATA cell:
% data = {yQ;nQ;yO;nO;nyes;nno;total_n};
% yQ = data{1};
% nQ = data{2};
% yO = data{3};
% nO = data{4};
% nyes = data{5};
% nno = data{6};
% N = data{7};
% extrt = data{8};

nneeded = 11;
%% Find deviance, by response
if data{6}<=nneeded, Dn = 0; % if number of zeros <12, ignore
else       Dn = fitanswer(param,data{4},data{2},0,data{7},data{8},method);
end

if data{5}<=nneeded, Dy = 0; % if number of ones <12, ignore
else       Dy = fitanswer(param,data{3},data{1},1,data{7},data{8},method);
end

D = Dy+Dn;

%% %%
%% Subfunction
function X2 = fitanswer(param,O,Q,reply,total_n,extrt,method)

%% Find expected proportions
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

if ~method
    % Calculate Chi-Square
    E = total_n*P;
    E(E<1e-5) = 1e-5;
    X2 = sum((O-E).^2./E);
else
    % Calculate -2*loglikelihood
    P(P<1e-5)=1e-5;
    X2 = -2*sum(log(P).*O);
end
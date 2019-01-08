%% FIT4PL  Fit a nonlinear model to a parameter (demo script)
% By default, DMAT only allows users to fit linear models (because it uses
% design matrices to apply simple models. However, this is only a practical
% restriction, and not a theoretical one. A simple feature (the
% |'NoFitting'| option of |multiestv4|) of DMAT can help more experienced
% users to fit models of higher complexity.

%% Generate some data with nonlinear drifts
% We're generating data where all parameters are constant, except for the
% drift rate, which varies nonlinearly with a covariate |X|. The nonlinear
% link function is the four-parameter logistic function:
%%
%
% $$v = l + (h-l) \frac{e^{aX+b}}{1+e^{aX+b}}$$
%
%%
% The parameters are:
%%
% * Upper asymptote _h_ = 0.30;
% * Lower asymptote _l_ = -0.15;
% * Location _b_ = 0.50;
% * Slope _a_ = 1.50;
clear all
close all
x = linspace(-5,4,10);
logistic = @(p,x) p(2)+(p(1)-p(2))./(1+exp(-p(4).*x-p(3)));
p = [0.3 -0.15 0.5 1.5]; % [h l b a]
drifts = logistic(p,x);
param = repmat(standardparset(0,0),5,1);
param(:,7) = drifts
%% 
% This is what the drifts look like:
plot(x,drifts,'k:')
axis square
box off
hold on
%%
data = multisimul(param,10000,1);
%% Build an OPTIONS structure as you usually would
% But add the field 'NoFitting' and assign it a true value (any integer
% except |0|, or just |true|).
desmat = [{'1'},{'1'},{'1'},{'1'},{'1'},{'1'},{[]}];
wd = 6; % number of free variables without drift, keeps code flexible
options = struct('DesignMatrix',{desmat},...
    'Name','Nonlinear Demo I',...
    'ShortSimplexRuns',5,...
    'MaxIter',10000,...
    'ObjectiveDecimals',4,...
    'NoFitting',true);
%% Use MULTIESTV4 to prepare the objective function
% Calling |multiestv4| with this options structure will not initiate
% parameter estimation, but will merely prepare the objective function and
% add it to the options structure:
options = multiestv4(data,options)
%% Adapt the objective to your needs
% Now, starting from this objective function, we can construct a new,
% higher-order objective function that reduces the number of parameters
% even more. For example:
newobj = @(y) options.objecfun([reshape(y(1:wd),1,wd),...
    logistic(y(wd+(1:4)),x)]);
%% Estimate parameters of the new model
% Now, rearrange the fields a bit, inserting the new objective:
options.oldobj = options.objecfun;
options.objecfun = newobj;
options.NoFitting = false;
%%
% Also insert a new initial guess (make sure you've got enough
% parameters!):
guess = [.5 -.1 0 1];
options.controls.small = [options.controls.small(1:wd)' guess];
options.controls.small([5 6])=[.03 .1]; % far away from the edge
%%
% And run the generative algorithm:
[output options]=genalg(options);

%% Process the output carefully
% The output as |genalg| returns it, is based on a linear model and thus
% not completely correct. You need to extract the design vector of the
% level 2 model and restore the level 1 design vector from that.
designvec_l2 = options.controls.small;
designvec_l1 = [designvec_l2(1:wd) logistic(designvec_l2(wd+(1:4)),x)];
options.controls.small = designvec_l1;
options.controls = bigger(options.controls);
minimum = options.controls.large;
%% Plot drift rates from the nonlinear model
plot(x,minimum(:,7),'rx-')
%% Estimate a free model
output2 = multiestv4(data,'DesignMatrix',desmat,...
    'Name','Nonlinear Demo II',...
    'ShortSimplexRuns',5,...
    'MaxIter',10000,...
    'ObjectiveDecimals',4,...
    'Guess',minimum);
%% Plot drift rates from the free model
plot(x,output2.Minimum(:,7),'bo')

%%
%%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.
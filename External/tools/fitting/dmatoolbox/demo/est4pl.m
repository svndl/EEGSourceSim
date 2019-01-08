%% EST4PL  Fit a nonlinear model to a parameter (demo script)
% By default, DMAT only allows users to fit linear models (because it uses
% design matrices to apply simple models. However, this is only a practical
% restriction, and not a theoretical one. With more advanced use of DMAT,
% models of higher complexity can be implemented. To do this, you need to
% bypass the |multiestv4| function as demonstrated in this script.

%% Generate some data with nonlinear drifts
% We're generating data where all parameters are constant, except for the
% drift rate, which varies nonlinearly with a covariate |X|. The nonlinear
% link function we use is the four-parameter logistic function:
%%
%
% $$v = l + (h-l) \frac{e^{bX+c}}{1+e^{bX+c}}$$
%
%%
% The parameters are:
%%
% * Upper asymptote _h_ = 0.30;
% * Lower asymptote _l_ = -0.15;
% * Location _c_ = 0.50;
% * Slope _b_ = 1.50;
%%
% First, clean up the workspace
clear all
close all
%%
% Generate the covariate |X|:
x = linspace(-5,4,10);
%%
% Prepare the four-parameter logistic link:
logistic = @(p,x) p(2)+(p(1)-p(2))./(1+exp(-p(4).*x-p(3)));
p = [0.30 -0.15 0.50 1.50]; % [h l c b]
%%
% Drift rates are a logistic function of the covariate:
drifts = logistic(p,x);
%%
% So the parameter set is this:
param = repmat(standardparset(0,0),5,1);
param(:,7) = drifts
%% 
% This is what the drifts look like:
plot(x,drifts,'k:')
axis square
box off
hold on
%%
% Now generate data from these parameters:
data = multisimul(param,10000,1);

%% First fit the free model
% There are two ways to implement nonlinear designs. The easy way,
% described in this demo, involves fitting the free model (i.e., no
% restrictions on drift rate across conditions) first, and then using the
% point estimates of the drift rates to get a good initial guess of the
% logistic parameters. But, there may be several reasons for which you
% might not want to fit this free model first, and then you should apply
% the method demonstrated in |fit4pl_hard|.
%%
% In this demo, we start by fitting the free model in the usual way:
desmat = [{'1'},{'1'},{'1'},{'1'},{'1'},{'1'},{[]}];
output.Minimum = param;
output = multiestv4(data,'DesignMatrix',desmat,...
    'Name','Free model');

%%
% Overlay the previous plot with recovered drift rates:
plot(x,output.Minimum(:,7),'bo')

%% Obtaining the link parameters
% Obtaining the four parameters of the link function is an optimization
% problem in itself. Fortunately, fitting a nonlinear regression line
% through the drift rates isn't very difficult. Here, we use the Statistics
% Toolbox function |nlinfit|. If you don't have it... Well, nonlinear
% regression isn't _that_ hard.
logpars = nlinfit(x,output.Minimum(:,7)',logistic,p)

%% Estimate parameters of the new model
% The previous execution of |multiestv4| has provided us with an objective
% function that is almost the one we need (in the field
% output.Options.objecfun), and with a design vector that minimizes this
% objective (in output.Options.controls.small). Starting from this
% objective function, we can construct a new, higher-order objective
% function that reduces the number of parameters even further. This
% second-level objective needs to accept ten parameters: one for each of
% the other diffusion model parameters, and four for the link function with
% which it will construct drift rates for all conditions.
%%
% We can formulate the objective like this, taking only one input vector,
% and transforming that into a vector that the old objective will like:
wd=6; % number of free parameters excluding drift, for reusability gcp
newobj = @(y) output.Options.objecfun([reshape(y(1:wd),1,wd),...
    logistic(y(wd+(1:4)),x)]);
%%
% Now we insert the new objective where the old one was, and store the old
% objective somewhere safe:
output.Options.oldobj = output.Options.objecfun;
output.Options.objecfun = newobj;
%%
% We obtained an initial guess for the link parameters and the other
% diffusion parameters from the free model. We put that initial guess into
% the output.Options.controls.small field:
guess = logpars(:);
output.Options.controls.small = [output.Options.controls.small(1:wd);guess];
%%
% And run the generative algorithm:
[ign options]=genalg(output.Options);

%% Process the output carefully
% The output as |genalg| returns it, is based on a linear model and thus
% not completely correct. We need to extract the design vector of the
% level 2 model and restore the level 1 design vector from that.
%%
% These are the link parameters (the seventh through tenth elements of the
% design vector):
linkpar = options.controls.small(wd+(1:4));
%%
% And this is how we recover the drift rates from the parameters:
designvec_l2 = options.controls.small;
designvec_l1 = [designvec_l2(1:wd);logistic(designvec_l2(wd+(1:4)),x)'];
options.controls.small = designvec_l1;
options.controls = bigger(options.controls);
minimum = options.controls.large
%%
% And overlay the recovered drifts in the plot again:
plot(x,minimum(:,7),'rx-')
hold off

%%
%%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.
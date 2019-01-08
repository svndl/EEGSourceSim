%% DMAT: Tips and tricks
% This is a list of running tips and cunning tricks that might help you get
% the most out of our software. Use them wisely, and they may serve you
% well.
%% Try to define models such that all parameters have similar ranges
% MATLAB's optimizers use the same step size and sensitivity for each
% parameter. That means that they are most effective if all parameters are
% on the same scale. In DMAT, the standard parameters have a range
% somewhere between 0 and 1. So, if you define designs, try to make sure
% that reasonable parameter estimates are also in that range.
%%
% For example, if you define a design matrix on the basis of a covariate,
% like this:
%%
%  >> E = 0:20:100
E = 0:20:100
%%
%  >> designmatrix = [ones(length(E),1) E']
designmatrix = [ones(length(E),1) E']
%%
% You should rescale it, such that most of the values are in the range 0-1,
% or maybe even lower (like 0-.5):
%%
%  >> E = E./max(E(:))
E = E./max(E(:))
%%
%  >> designmatrix = [ones(length(E),1) E']
designmatrix = [ones(length(E),1) E']
%% Try different starting points
% In particular if you have more complicated designs, the EZDIFF starting
% point may not be ideal (though it usually is very good). Restarting the
% optimization procedure with a different initial guess might yield better
% results.

%% Author of this file
% Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
% Part of the DMA Toolbox. Please read the End User License Agreement,
% contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
% See also http://ppw.kuleuven.be/okp/dmatoolbox.
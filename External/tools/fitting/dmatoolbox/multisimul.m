function data = multisimul(par,N,seed)
%MULTISIMUL  Generates data for multiple conditions from a diffusion model
%   DATA = MULTISIMUL(PARAM,N), where PARAM is a valid parameter set, N is
%   a vector carrying the number of data points in each condition (or a
%   scalar, in which case all conditions will have the same number of data
%   points), returns DATA, an N-by-3 data matrix with condition in the
%   first column, binary responses in the second and reaction times (in
%   seconds) in the third column.
%
%   DATA = MULTISIMUL(PARAM,N,SEED), uses SEED to start the random
%   generator.
%
%   DATA = MULTISIMUL(PARAM), uses all N = 250.
%
%   DATA = MULTISIMUL, uses the standard parameter set 1, as contained in
%   STANDARDPARSET.
%
%   The volatility parameter s is fixed to .1.
%
%   See also SIMULDIFF, STANDARDPARSET.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Added optional use of MEX-file SIMDIF (not included in DMAT).

if nargin<3
    seed = [];
    if nargin<2
        N=250;
        if nargin<1
            par = standardparset(1);
        end
    end
end

r = size(par,1);
if length(N)==1
    N = repmat(N,r,1);
elseif length(N)~=r
    error('Number of parameter sets does not agree with size of N.');
end

data = [];

if ~isempty(seed)
    rand('seed',seed);
    randn('seed',seed);
end

for ctr=1:r
    try
        [T,XX] = simdif(par(ctr,:),N(ctr),ceil(rand*intmax));
    catch
        [T,XX] = simuldiff(par(ctr,:),N(ctr));
    end
    data = [data;repmat(ctr,N(ctr),1),XX',T'];
end

% This shuffles the trials, no need for it beyond testing DMAT
[ignore,p] = sort(rand(1,length(data))); %#ok
data = data(p,:);
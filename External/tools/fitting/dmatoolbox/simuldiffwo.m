function [T,XX,fastOlInd,slowOlInd] = simuldiffwo(par,N,fastolprop,...
    fastolrange,slowolprop,slowolrange)
% SIMULDIFFWO  Generates data according to a diffusion model, with outliers
%   [T,XX,FASTOLIND,SLOWOLIND] = SIMULDIFFWO(PARAM,N,FASTOLPROP,...
%       FASTOLRANGE,SLOWOLPROP,SLOWOLRANGE), where PAR and N are inputs to
%   SIMULDIFF (see help for that function), and FASTOLPROP and SLOWOLPROP
%   define the proportions of fast and slow outliers, respectively, and
%   FASTOLRANGE and SLOWOLRANGE define their range (they follow a uniform
%   distribution over this range), returns T, a 1-by-N vector with reaction
%   times (in seconds), and X, a 1-by-N vector with binary responses. These
%   data will contain RT outliers. Fast outliers are guesses (accuracy
%   approximately .50) and slow outliers are 'delayed startups' (accuracy
%   approximately equal to normal observations).
%
%   See also MULTISIMUL, SIMULDIFF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<6
    error('DMAT:simuldiffwo:notEnoughInputs',...
        'SIMULDIFFWO requires exactly six parameters.')
end

[T XX] = simuldiff(par,N);
slowOlN = round(slowolprop*N);
fastOlN = round(fastolprop*N);

olN = slowOlN+fastOlN;
olInd = randperm(olN);
slowOlInd = olInd(1:slowOlN);
fastOlInd = olInd(slowOlN+1:olN);

slowOls = rand(1,slowOlN)*(slowolrange(2)-slowolrange(1))+slowolrange(1);
fastOls = rand(1,fastOlN)*(fastolrange(2)-fastolrange(1))+fastolrange(1);

T(slowOlInd)=slowOls;
T(fastOlInd)=fastOls;
XX(fastOlInd)=round(rand(size(fastOlInd)));

if ~nargout
    fprintf(1,['\n  Simulation summary:\n\n  Mean(X) = %f\n  RT(X=1) ',...
        '= %f\n  RT(X=0) = %f\n\n'],mean(XX),mean(T(XX==1)),mean(T(XX==0)))
end
if nargout == 1
    T = [XX' T'];
end
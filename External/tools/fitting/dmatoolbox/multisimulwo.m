function [data,folind,solind] = multisimulwo(par,N,seed,fastolprop,...
    fastolrange,slowolprop,slowolrange)
%MULTISIMULWO  Sample multiple conditions w/ outliers from diffusion model
%   [DATA,FOLIND,SOLIND]= MULTISIMULWO(PARAM,N,SEED,FASTOLPROP,...
%      FASTOLRANGE,SLOWOLPROP,SLOWOLRANGE). See MULTISIMUL and SIMULDIFFWO
%   for usage. 
%
%   See also MULTISIMUL, SIMULDIFFWO, and STANDARDPARSET.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

r = size(par,1);
if length(N)==1;
    N = N*ones(r,1);
elseif length(N)~=r;
    error('Number of parameter sets does not agree with size of N.');
end

data = [];

if nargin>=3 && ~isempty(seed);
    rand('seed',seed);
    randn('seed',seed);
end

folind=[];
solind=[];
for ctr=1:r
    if nargin>3 && (fastolprop || slowolprop),
        [T,XX,folindn,solindn] = simuldiffwo(par(ctr,:),N(ctr),fastolprop,fastolrange,slowolprop,slowolrange);
    else
        [T,XX] = simuldiff(par(ctr,:),N(ctr));
        folindn=[];
        solindn=[];
    end
    folind=[folind(:);folindn(:)+(ctr-1)*length(T)];
    solind=[solind(:);solindn(:)+(ctr-1)*length(T)]; 
    data = [data;repmat(ctr,N(ctr),1),XX',T'];
end
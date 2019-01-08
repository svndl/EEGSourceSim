function o = pctile(d,p)
% PCTILE  Compute percentiles of a vector
%   O = PCTILE(D,P), where D is a vector of data and P a vector of
%   percentiles, returns O, a vector of associated quantiles.
%
%   See also: PRCTILE (Statistics Toolbox).
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Added in 0.4 to remove dependence on Statistics Toolbox.

try
    o = prctile(d,p);
catch
    n = length(d);
    d = sort(d);
    spc = [0 100*(0.5:(n-0.5))./n 100]';
    vec = [d(1,:); d(1:n,:); d(n,:)];
    o = interp1q(spc,vec,p')';
end
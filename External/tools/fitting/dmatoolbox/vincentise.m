function vinc = vincentise(RTcell,pct)
% VINCENTISE  Performs quantile averaging of reaction time data
%   CDF = VINCENTISE(RTCELL,PCT), where RTCELL is a cell array containing
%   one participant's RTs in each cell. PCT is a vector containing
%   percentiles at which to center. 
%
%   CDF = VINCENTISE(RTCELL), centers at .05:.05:.95.
%
%   (Note: This procedure is not actually Vincentization, but the name is
%   kept for backward compatibility reasons. Future versions of DMAT will
%   likely use the quantile averaging procedure in a more integrated
%   fashion.)
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: No longer crashes with too few elements in a data set.

%% Preparation
if nargin<2, pct=.05:.05:.95; end

np = numel(RTcell);
tcdf = zeros(length(pct),np);
keeps = true(np,1);
nneeded = 11;

for ctr=1:np
    t = sort(RTcell{ctr});
    if length(t)>nneeded
        tcdf(:,ctr) = t(ceil(pct*length(t)));
    else
        keeps(ctr) = false;
    end
end

%% Vincentizing
vinc = mean(tcdf(:,keeps), 2);

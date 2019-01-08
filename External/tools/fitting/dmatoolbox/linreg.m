function b = linreg(y,x)
% LINREG  Linear regression of X on Y
%   B = LINREG(Y,X), where Y is a set of diffusion parameters and X is a
%   design matrix, returns B, a least-squares estimate of the model
%   parameter vector.
%
%   See also: REGRESS (Statistics Toolbox).
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Added in 0.3 to remove dependence on Statistics Toolbox.

try
    b = regress(y,x);
catch
    ms = 'Statistics Toolbox not found. DMAT may work at lower precision.';
    lw = lastwarn;
    if ~strcmp(ms,lw)
        warning('DMAT:linreg:StatisticsToolboxNotFound',ms)
    end
    b = pinv(x)*y;
end
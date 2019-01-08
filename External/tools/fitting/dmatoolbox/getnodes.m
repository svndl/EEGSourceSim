function varargout = getnodes
%GETNODES  Tells you how many nodes are currently used by CDFDIF
%   GETNODES() prints the number of nodes used by CFDIF to the screen.
%
%   [NR_V NR_Z] = GETNODES() returns the number of quadrature nodes for the
%   drift rate integral (Gauss-Hermite nodes) to NR_V, and the number of
%   nodes for the starting point integral (Gauss-Legendre nodes) to NR_Z.
%
%   See also SETNODES, CDFDIF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

nr = getpref('dmatoolbox','nodes');
if ~nargout
    fprintf(['Binary file is cdfdif.%s.\nGauss-Hermite quadrature '...
        'points: %3i\nGauss-Legendre quadrature points:%3i\n'],...
        mexext,nr(1),nr(2));
else
    varargout = num2cell(nr);
end
function p=isposdef(h,tol)
%ISPOSDEF  Check if input matrix is positive definite
%   P = ISPOSDEF(H), returns 1 if all eigenvalues of H are larger than zero
%   and 0 otherwise.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<2, tol=0; end
p = ndims(h)==2 && ~diff(size(h)) && ~any(isnan(h(:))) && all(eig(h)>-tol);
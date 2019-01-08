function b=isnested(A,B)
%ISNESTED  Check if A is nested in B.
%   BOOL = ISNESTED(A,B), returns 1 (true) if design matrix A is nested in
%   design matrix B and 0 (false) otherwise.
%
%   WARNING: THIS FUNCTION MAY RETURN FALSE POSITIVES. DO NOT RELY ON IT.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

% TODO % Something something like this might work better:
% b = isequal(orth(B)*orth(B')*A,A)

[r c]=size(A);
[s d]=size(B);
b = false;

% If NaNs, return NaN and bail out
if any([isnan(A(:));isnan(B(:))])
    b = NaN;
    return
end

% If A or B are not of full rank, return false
if rank(A)~=c || rank(B)~=d
    return
end

% If B smaller than A, A can't be nested in it
if ~(s<r || d<c)
    C=rref(A);
    D=rref(B);
    b=isequal(D(1:r,1:c),C);
end

function b = iseye(a)
%ISEYE  Check if input is an identity matrix
%   B = iseye(A), returns 1 (true) if A is an identity matrix and 0 (false)
%   otherwise.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

[r c] = size(a);
b = isnumeric(a) && (r==c) && isequal(a,eye(r)) && ~any(isnan(a(:)));

% Might be faster?
% b = ~any(diff(size(a)))&&all(sum(a)==sum(a')==diag(a)')
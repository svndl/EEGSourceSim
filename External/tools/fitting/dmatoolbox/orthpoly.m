function c=orthpoly(cond,degree)
%ORTHPOLY  Table of coefficients of orthogonal polynomials
%   COEFF = ORTHPOLY(C,D) returns coefficients for a polynomial of degree D
%   (1-5), suitable for a design with C (3-10) conditions. Can be useful
%   for constructing design matrices. COEFF is a column matrix.
%   If D is a vector, a design matrix with those orders is returned.
%
%   Note that if you use these polynomials in a design matrix, it is
%   advisable to rescale them so that the basic parameters are in the range
%   [0 .5] or [0 1].
%
%   Not all combinations of C and D are available, but the m-file is easy
%   to expand.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

% TODO % There's an algorithm for doing these, hidden in SETNODES
% somewhere.

%% Initialize output
c = zeros(cond,length(degree));

%% If no conditions, return empty design matrix
if cond<=0
    return
end

%% If oversaturated, throw error
if (~isscalar(degree) && any(degree>=cond)) || ...
        (isscalar(degree) && degree>cond)
    error('DMAT:orthpoly:youAskTooMuch','This design is oversaturated.')
end

%% If multiple degrees were requested, cycle through them
if ~isscalar(degree)
    for ctr = 1:length(degree)
        c(:,ctr) = orthpoly(cond,degree(ctr));
    end
    return
end

%% Zeroeth degree is a column of ones
if degree==0
    c = ones(cond,1);
    return
end

%% Saturated designs are always equivalent to identity matrices
if cond==degree
    c = eye(cond);
    return
end

%% Anything else, check the table
% This part could be expanded with more orthogonal polynomials
switch cond
    case 3
        switch degree
            case 1, c = [-1 0 1];                
            case 2, c = [1 -2 1];
        end
    case 4
        switch degree
            case 1, c = [-3 -1 1 3];
            case 2, c = [1 -1 -1 1];
            case 3, c = [-1 3 -3 1];
        end
    case 5
        switch degree
            case 1, c = [-2 -1 0 1 2];
            case 2, c = [2 -1 -2 -1 2];
            case 3, c = [-1 2 0 -2 1];
            case 4, c = [1 -4 6 -4 1];
        end
    case 6
        switch degree
            case 1, c = [-5 -3 -1 1 3 5];
            case 2, c = [5 -1 -4 -4 -1 5];
            case 3, c = [-5 7 4 -4 -7 5];
            case 4, c = [1 -3 2 2 -3 1];
        end
    case 7
        switch degree
            case 1, c = [-3 -2 -1 0 1 2 3];
            case 2, c = [5 0 -3 -4 -3 0 5];
            case 3, c = [-1 1 1 0 -1 -1 1];
            case 4, c = [3 -7 1 6 1 -7 3];
        end
    case 8
        switch degree
            case 1, c = [-7 -5 -3 -1 1 3 5 7];
            case 2, c = [7 1 -3 -5 -5 -3 1 7];
            case 3, c = [-7 5 7 3 -3 -7 -5 7];
            case 4, c = [7 -13 -3 9 9 -3 -13 7];
            case 5, c = [-7 23 -17 -15 15 17 -23 7];
        end
    case 9
        switch degree
            case 1, c = [-4 -3 -2 -1 0 1 2 3 4];
            case 2, c = [28 7 -8 -17 -20 -17 -8 7 28];
            case 3, c = [-14 7 13 9 0 -9 -13 -7 14];
            case 4, c = [14 -21 -11 9 18 9 -11 -21 14];
            case 5, c = [-4 11 -1 -9 0 9 4 -11 4];
        end
    case 10
        switch degree
            case 1, c = [-9 -7 -5 -3 -1 1 3 5 7 9];
            case 2, c = [6 2 -1 -3 -4 -4 -3 -1 2 6];
            case 3, c = [-42 14 35 31 12 -12 -31 -355 -14 42];
            case 4, c = [18 -22 -17 3 18 18 3 -17 -22 18];
            case 5, c = [-6 14 -1 -11 -6 6 11 1 -14 6];
        end
end
%% End table
% Do not edit beyond this

%% If combination not in table, c is still filled with zeros
if ~any(c)
    error('DMAT:orthpoly:CombinationNotInTable',...
        'The combination C=%i, D=%i is not in the table.',cond,degree)
else
    c=c(:);
end
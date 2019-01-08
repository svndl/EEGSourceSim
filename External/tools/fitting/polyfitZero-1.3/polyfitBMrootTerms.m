function [p,S,mu] = polyfitBMrootTerms(x,y,degree,b,m,root,terms)
% POLYFITBM Fit polynomial to data, forcing y-intercept, slope and root
%   P = POLYFITBM(X,Y,N,B,M,ROOT,TERMS) is similar POLYFIT(X,Y,N) except that the
%   y-intercept is forced to B, i.e. P(N) = B, and the slope at X=0 is forced
%   to M, i.e. DY(X=0)/DX=M. In the same way as POLYFIT, the coefficients,
%   P(1:N-1), fit the data Y best in the least-squares sense. You can also use
%   Y = POLYVAL(PZERO,X) to evaluate the polynomial because the output is the
%   same as POLYFIT.
%
%   [P,S,MU] = POLYFITBM() Return structure, S, similar to POLYFIT for use
%   with POLYVAL to calculate error estimates of predictions with P.
%
%   [P,S,MU] = POLYFITBM() Scale X by std(X), returns MU = [0,std(X)].
%
%   See also POLYVAL, POLYFIT
%
%   Also see <a href="http://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn">POLYFITN by John D'Errico</a>
%

% Copyright (c) 2013 Mark Mikofski
% Version 1-1, 2013-10-15
%   add delta output
%   center and scale
% Version 1-0, 2011-06-29

%% check args
% X & Y should be numbers
assert(isnumeric(x) && isnumeric(y),'polyfitZero:notNumeric', ...
    'X and Y must be numeric.')
dim = numel(x); % number of elements in X
% DEGREE should be scalar positive number between 1 & 10 inclusive
assert(isnumeric(degree) && isscalar(degree) && degree>0 && degree<=10, ...
    'polyfitZero:degreeOutOfRange', ...
    'DEGREE must be an integer between 1 and 10.')
% DEGREE must be less than number of elements in X & Y
assert(degree<dim && degree==round(degree), ...
    'polyfitZero:DegreeGreaterThanDim', 'DEGREE must be less than numel(X)')
% X & Y should be same size vectors
assert(isvector(x) && isvector(y) && dim==numel(y), ...
    'polyfitZero:vectorMismatch', 'X and Y must be vectors of the same length.')
assert(isnumeric(b) && isscalar(b),'polyfitZero:badIntercept', ...
    'B, y-intercept, must be scalar number.')
assert(isnumeric(m) && isscalar(m),'polyfitZero:badSlope', ...
    'M, slope, must be scalar number.')
if nargin<4 || isempty(b), b=0;end % set default intercept to zero, force to origin
if nargin<5 || isempty(m), m=0;end % set default slope to zero, force no slope
%% solve
% convert X & Y to column vectors
x = x(:); y = y(:);
% Scale X.
% attribution: this is based on code from POLYFIT by The MathWorks Inc.
if nargout > 2
   mu = [0; std(x)];
   x = (x - mu(1))/mu(2);
end
% using pow() is actually as fast or faster than looping, same # of flops!
nterms = numel(terms);
z = zeros(dim,numel(nterms));
idx = 1;
for n = terms
    z(:,idx) = x.^(degree-n);
    idx = idx+1;
end
b_root = -b/root;
m_root = -(m-b_root)/root;
y_root = y./(x-root) - b_root - m_root.*x;
p_root = z\y_root; % solve
p_terms = zeros(nterms,1);
p_terms(terms) = p_root;
p_terms = [p_terms;m_root;b_root];
p = conv(p_terms,[1,-root]);
%% error estimates
% attribution: this is based on code from POLYFIT by The MathWorks Inc.
if nargout > 1
    V = [z,x,ones(dim,1)]; % append constant term for Vandermonde matrix
    % Return upper-triangular factor of QR-decomposition for error estimates
    R = triu(qr(V,0));
    r = y - V*p;
    S.R = R(1:size(R,2),:);
    S.df = max(0,length(y) - (degree+1));
    S.normr = norm(r);
end
p = p'; % polynomial output is row vector by convention
end
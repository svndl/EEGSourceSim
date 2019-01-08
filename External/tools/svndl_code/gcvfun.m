function G = gcvfun(lambda,s2,beta,delta0,mn,t_nb)

% Auxiliary routine for gcv.  PCH, IMM, 12/29/97.

if (nargin==6)
   f = (lambda^2)./(s2 + lambda^2);
else
   f = lambda./(s2 + lambda);
end

f_new = repmat(f,t_nb,1);
G = (norm(f_new.*beta)^2 + delta0)/ ((mn + sum(f))^2);
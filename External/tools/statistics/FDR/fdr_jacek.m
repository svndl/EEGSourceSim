function h = fdr(pvals,q)
% pvals: p-values of your multiple significance tests
% q: desired FDR level (like alpha) -- defaults to 0.05
% h: array of size(pvals) indicating the corrected hypothesis test results:
% 1 indicates a rejected null hypothesis
% reference: Benjamini, Yoav; Yekutieli, Daniel (2001). "The control of the false discovery rate in multiple testing under dependency". Annals of Statistics 29 (4): 1165–1188.
sz=size(pvals); pvals=pvals(:);  m=length(pvals);
if size(nargin)<2, q=0.05; end
[c,d]=sort(pvals,'ascend');
mxind=max(find( (c<=(1:m)'/m*q) ));
h=zeros(size(pvals));
h(d(1:mxind))=1;  
h=reshape(h,sz(1),sz(2));



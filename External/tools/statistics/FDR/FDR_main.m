% General FDR procedure
% Parameters:
%   P         p-value vector
%   q         scalar controlling FDR output
%   str    name of procedure (case insensetive)
%   Output:
%   F        vector containing indices of hypotheses that passsed selected procedure 

function F=FDR_main(P, q, str)
  switch lower(str)
   case 'ibh_up'
      F = fdr_IBH_step_up(P,q);
   case 'ibh_down'
      F = fdr_IBH_step_down(P,q);
   case 'bh95'
      F = fdr_proc(P,q);
   case 'sth'
       F = fdr_STH(P,q);
    case 'bky'
        F = fdr_BKY(P,q);
   otherwise
      disp('Unknown procedure.')
  end
end


%FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a p-values vector (up to 30000 elements long), q is a scalar which
%control the false discovery rate. output: F vector contain the index of
%the hypotheses that pass the procedure. Improve Benjamini Hochberg 1995
function F=fdr_IBH_step_down(P,q)

[n,m]=size(P);
if n>m; P=P';m=n;end
psor=sort(P);
m_C_s=1e5 *[     0.001000000000000   0.000010306035633   0.000350000000000
   0.002000000000000   0.000010199148511   0.000550000000000
   0.003000000000000   0.000010156707816   0.000720000000000
   0.004000000000000   0.000010132670482   0.000860000000000
   0.005000000000000   0.000010117091093   0.000980000000000
   0.006000000000000   0.000010105538438   0.001090000000000
   0.007000000000000   0.000010096875092   0.001190000000000
   0.008000000000000   0.000010090004408   0.001290000000000
   0.009000000000000   0.000010084411050   0.001380000000000
   0.010000000000000   0.000010079681966   0.001470000000000
   0.020000000000000   0.000010054904854   0.002170000000000
   0.030000000000000   0.000010044257282   0.002720000000000
   0.040000000000000   0.000010038076350   0.003180000000000
   0.050000000000000   0.000010033856638   0.003590000000000
   0.060000000000000   0.000010030847440   0.003960000000000
   0.070000000000000   0.000010028436400   0.004300000000000
   0.080000000000000   0.000010026535868   0.004620000000000
   0.090000000000000   0.000010025022004   0.004910000000000
   0.100000000000000   0.000010023656804   0.005210000000000
   0.150000000000000   0.000010019218320   0.006450000000000
   0.200000000000000   0.000010016616056   0.007500000000000
   0.250000000000000   0.000010014818536   0.008430000000000
   0.300000000000000   0.000010013492425   0.009280000000000
   0.400000000000000   0.000010011679484   0.010770000000000
   0.500000000000000   0.000010010408058   0.012110000000000
   0.600000000000000   0.000010009486880   0.013320000000000
   0.700000000000000   0.000010008794243   0.014390000000000
   0.800000000000000   0.000010008212513   0.015430000000000
   0.900000000000000   0.000010007735493   0.016410000000000
   1.000000000000000   0.000010007339726   0.017310000000000];
C=interp1(m_C_s(:,1),m_C_s(:,2),m,'spline');
s=interp1(m_C_s(:,1),m_C_s(:,3),m,'spline');
m0_hat=C*min(m,max(s,2*sum(P)));%the estimator
F=zeros(n,1);
fdr_line=[1:m]*q/m;


in_fdr=find((psor-fdr_line*m/m0_hat)>0,1,'first');
if in_fdr>1
    F=find(P<=psor(in_fdr-1));    
elseif isempty(in_fdr)==1
    F=[1:m];
else
    F=[];
end
end

%FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a matrix of sortdet p-values (each row is a vector), q is a scalar which
%control the false discovery rate. output: F vector contain the number of
%hypotheses that pass the procedure. Improve Benjamini Hochberg 1995
function F=fdr_IBH_step_up(P,q)

[n,m]=size(P);
if n>m; P=P';m=n;end
psor=sort(P);
m_C_s=1e5 *[     0.001000000000000   0.000010306035633   0.000350000000000
   0.002000000000000   0.000010199148511   0.000550000000000
   0.003000000000000   0.000010156707816   0.000720000000000
   0.004000000000000   0.000010132670482   0.000860000000000
   0.005000000000000   0.000010117091093   0.000980000000000
   0.006000000000000   0.000010105538438   0.001090000000000
   0.007000000000000   0.000010096875092   0.001190000000000
   0.008000000000000   0.000010090004408   0.001290000000000
   0.009000000000000   0.000010084411050   0.001380000000000
   0.010000000000000   0.000010079681966   0.001470000000000
   0.020000000000000   0.000010054904854   0.002170000000000
   0.030000000000000   0.000010044257282   0.002720000000000
   0.040000000000000   0.000010038076350   0.003180000000000
   0.050000000000000   0.000010033856638   0.003590000000000
   0.060000000000000   0.000010030847440   0.003960000000000
   0.070000000000000   0.000010028436400   0.004300000000000
   0.080000000000000   0.000010026535868   0.004620000000000
   0.090000000000000   0.000010025022004   0.004910000000000
   0.100000000000000   0.000010023656804   0.005210000000000
   0.150000000000000   0.000010019218320   0.006450000000000
   0.200000000000000   0.000010016616056   0.007500000000000
   0.250000000000000   0.000010014818536   0.008430000000000
   0.300000000000000   0.000010013492425   0.009280000000000
   0.400000000000000   0.000010011679484   0.010770000000000
   0.500000000000000   0.000010010408058   0.012110000000000
   0.600000000000000   0.000010009486880   0.013320000000000
   0.700000000000000   0.000010008794243   0.014390000000000
   0.800000000000000   0.000010008212513   0.015430000000000
   0.900000000000000   0.000010007735493   0.016410000000000
   1.000000000000000   0.000010007339726   0.017310000000000];
C=interp1(m_C_s(:,1),m_C_s(:,2),m,'spline');
s=interp1(m_C_s(:,1),m_C_s(:,3),m,'spline');
m0_hat=C*min(m,max(s,2*sum(P)));%the estimator
F=zeros(n,1);
fdr_line=[1:m]*q/m;


in_fdr=find((psor-fdr_line*m/m0_hat)<=0,1,'last');
if in_fdr>0
    F=find(P<=psor(in_fdr));    
elseif isempty(in_fdr)==1
    F=[];
end
end

%FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a vector of p-values, q is a scalar which
%control the false discovery rate. output: F the indices of the element
%that pass the procedure
function F=fdr_proc(P,q)

[m,n]=size(P);
if n>m; P=P'; end
N=length(P);
[P_sor]=sort(P);
in_fdr=find((P_sor-[1:N]'*q/N)<=0,1,'last');
if in_fdr>0
    p_d=P_sor(in_fdr);
    F=find(P<=p_d);
elseif isempty(in_fdr)==1
    F=[];
else
    F=[1:N];
end
end

%FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a vector of p-values, q is a scalar which
%control the false discovery rate. output: F the indices of the element
%that pass the procedure. Storey, Taylor, Siegmund 2004
function F=fdr_STH(P,q)

[m,n]=size(P);
if n>m; P=P'; end
N=length(P);

m0=(N+1-length(find(P<0.5)))/0.5;

[P_sor]=sort(P);
% % %step down from the smaller p-value and up
in_fdr=find((P_sor-[1:N]'*q/m0)>0,1,'first');
if in_fdr>1
    p_d=P_sor(in_fdr-1);
    F=find(P<=p_d);
elseif isempty(in_fdr)==1
    F=[1:N];
else
    F=[];
end
end
%FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a vector of p-values, q is a scalar which
%control the false discovery rate. output: F the indices of the element
%that pass the procedure. Benjamini Krieger Yekutieli 2006
function F=fdr_BKY(P,q)

[m,n]=size(P);
if n>m; P=P'; end
[P_sor]=sort(P);
N=length(P);
in_fdr=find((P_sor-[1:N]'.*q./(N+1-[1:N]'*(1-q)))>0,1,'first');
if in_fdr>1
    p_d=P_sor(in_fdr-1);
    F=find(P<=p_d);
elseif isempty(in_fdr)==1
    F=[1:N];
else
    F=[];
end
end

%% polyfitB example
%% LaTex
%
% $$f\left(x\right)=p_1 x^n + p_2 x^{\left(n-1\right)} + \ldots + p_n x + b$$
%
%% initialize workspace
close('all'),clear('all'),clc
fprintf('polyfitB example: fit polynomial forcing intercept to B\n\n')
%% create some data with noise
npts = 11;xmax = 10;
intercept = 1+rand(1,1)/10;
fprintf('actual intercept: %g\n', intercept)
x = linspace(0,xmax,npts);y = intercept+(x+rand(1,npts)/10).^2;
%% fit data
degree = 2;b = 1;
fprintf('fit to degree: %g & intercept: %g\n', degree, b)
p = polyfitB(x,y,degree,b);
for n = 1:degree+1,fprintf('p%d = %f\n',n,p(n)),end
%% scale data
[p,~,mu] = polyfitB(x,y,degree,b);
fprintf('\nScale X:\n')
for n = 1:degree+1,fprintf('p%d = %f\n',n,p(n)),end
fprintf('scaled by %f\n',mu(2))
%% get error estimates
[p,S,mu] = polyfitB(x,y,degree,b);
[yest,derr] = polyval(p,x,S,mu); % fit to data, calculate error
plot(x,y,'o'),hold('all'),grid
errorbar(x,yest,derr),title('Polynomial fit forcing y through b.')
xlabel('x'),ylabel('y'),legend('data','fit','Location','NorthWest')
%% annotate
pos = get(gca,'Position');
xl = ceil(max(x)/(xmax/2))*(xmax/2);xlim([0,xl])
yl = ceil(max(y)/(xmax/2))*(xmax/2);ylim([0,yl])
xtrim = -0.05;
for n = 1:numel(x)
    xpos = min(max(0,pos(1)+pos(3)*x(n)/xl+xtrim),1);
    ypos = min(max(0,pos(2)+pos(4)*yest(n)/yl),1);
    annotation('textbox',[xpos,ypos,0.1,0.1], ...
        'LineStyle','none','FontWeight','bold', ...
        'String',sprintf('%4.2f%%',derr(n)/yest(n)*100))
end
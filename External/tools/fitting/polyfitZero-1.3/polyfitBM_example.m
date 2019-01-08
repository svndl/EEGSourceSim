%% polyfitBM example
%% LaTex
%
% $$f\left(x\right)=p_1 x^n + p_2 x^{\left(n-1\right)} + \ldots + m x + 1$$
%
%% initialize workspace
close('all'),clear('all'),clc
fprintf('polyfitB0 example: fit polynomial forcing intercept to B and slope to M\n\n')
%% create some data with noise
npts = 11;xmax = 10;
intercept = 1+rand(1,1)/10;
x = linspace(0,xmax,npts);y = (x+rand(1,npts)/10).^2;
theta = rand(1,1)+15; % [degrees] random rotation
fprintf('actual slope = %g, actual intercept = %g\n',tand(theta),intercept);
xyrot = [cosd(theta),sind(theta);-sind(theta),cosd(theta)]*[x;y]; % rotate randomly
x = xyrot(1,:);y = intercept+xyrot(2,:);
%% fit data
degree = 2;b = 1;m = tand(15);
fprintf('fit to degree: %g, intercept: %g & slope: %g\n', degree, b, m)
p = polyfitBM(x,y,degree,b,m);
for n = 1:degree+1,fprintf('p%d = %f\n',n,p(n)),end
%% scale data
[p,~,mu] = polyfitBM(x,y,degree,b,m);
fprintf('\nScale X:\n')
for n = 1:degree+1,fprintf('p%d = %f\n',n,p(n)),end
fprintf('scaled by %f\n',mu(2))
%% get error estimates
[p,S,mu] = polyfitBM(x,y,degree,b,m);
[yest,derr] = polyval(p,x,S,mu); % fit to data, calculate error
plot(x,y,'o'),hold('all'),grid
errorbar(x,yest,derr),title('Polynomial fit forcing y through b with slope m.')
xlabel('x'),ylabel('y'),legend('data','fit','Location','NorthWest')
%% annotate
pos = get(gca,'Position');
xl = [floor(min(x)/(xmax/2))*(xmax/2),ceil(max(x)/(xmax/2))*(xmax/2)];xlim(xl)
yl = [floor(min([y,yest])/(xmax/2))*(xmax/2),ceil(max([y,yest])/(xmax/2))*(xmax/2)];ylim(yl)
xtrim = -0.0;
for n = 1:numel(x)
    xpos = min(max(0,pos(1)+pos(3)*abs((x(n)-xl(1))/diff(xl))+xtrim),1);
    ypos = min(max(0,pos(2)+pos(4)*abs((yest(n)-yl(1))/diff(yl))),1);
    annotation('textbox',[xpos,ypos,0.1,0.1], ...
        'LineStyle','none','FontWeight','bold', ...
        'String',sprintf('%4.2f%%',derr(n)/yest(n)*100))
end
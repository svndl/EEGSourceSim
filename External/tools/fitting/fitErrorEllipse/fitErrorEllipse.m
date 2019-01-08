function [ampDiff,phaseDiff,zSNR,errorEllipse] = fitErrorEllipse(xyData,ellipseType,makePlot,returnRad)
%   [ampDiff,phaseDiff,zSNR,errorEllipse] = fitErrorEllipse(xyData,[ellipseType],[makePlot])
%   user provides xyData, an Nx2 matrix of 2D data of N samples
%   opt. input ellipseType can be 'SEM' '95CI' or a string specifying
%       a different percentage CI formated following: '%.1fCI'. Default is
%       'SEM'.
%   opt. input makePlot is a logical specifying whether or not to
%       generate a plot of the data & ellipse & eigen vectors (draw at
%       length of 1 std dev, i.e. sqrt(corresponding eigenvalue)). Default
%       is false.
%   opt. input returnRad is a logical specifying whether to return values
%       in radians or degrees. Default is false (degrees).
%
%   *NOTE: For EEG measurements in the Fourier domain, xyData rows should
%   be: [real,imag].
%
%   The function uses the eigen value decomposition to find the two
%   perpandicular axes aligned with the data in xyData, where the
%   eigenvalues of each eigen vector correspond to the variances along each
%   direction. An ellipse is then fit to this data, at a distance from the
%   mean datapoint depending on the type of ellipseType specified.
%
%   Calculations for the error ellipses based on alpha-specified confidence
%   regions (e.g. 95%CI or 68%CI) are calculated following information from
%   Chapter 5 of Johnson & Wickern (2007) Applied Multivariate Statistical
%   Analysis, Pearson Prentice Hall.
%
%   ampDiff and phaseDiff are the differences of the lower and uper bounds
%   from the mean amplitude and mean phase (mean - lower, upper - mean).
%
%   DEPENDENCY: eigFourierCoefs.m which is part of the GitHub repository
%   mrC in tools/fitting/fitErrorEllipse

xyData = double(xyData); 

if nargin<2 || isempty(ellipseType)
    ellipseType = 'SEM';
end
if nargin<3
    makePlot = false;
end
if nargin<4
    returnRad = false;
end

dims = size(xyData);
N = dims(1);
if dims(2) ~= 2
    error('input data must be a matrix of 2D row samples');
end
if N < 2
    error('input data must contain at least 2 samples');
end

srIx = 1;
siIx = 2;

try
    [meanXy,~,smaller_eigenvec,smaller_eigenval,larger_eigenvec,larger_eigenval,phi] = eigFourierCoefs(xyData);
catch
    fprintf('The eigen value decomposition of xyData could not be run, probably your data do not contain >1 sample.');
end

theta_grid = linspace(0,2*pi);

switch ellipseType
    case '1STD'
        a = sqrt(larger_eigenval); 
        b = sqrt(smaller_eigenval);
    case '2STD'
        a = 2*sqrt(larger_eigenval); 
        b = 2*sqrt(smaller_eigenval);
    case 'SEMarea'
        a = sqrt(larger_eigenval/sqrt(N)); % scales the ellipse's area by sqrt(N)
        b = sqrt(smaller_eigenval/sqrt(N));
    case {'SEMellipse' 'SEM'} % default
        a = sqrt(larger_eigenval)/sqrt(N); % contour at stdDev/sqrt(N)
        b = sqrt(smaller_eigenval)/sqrt(N);
    case '95CI'
        % following Eqn. 5-19 of Johnson & Wichern (2007):
        t0_sqrd = ( (N-1)*2 ) / ( N*(N-2) ) * finv( 0.95, 2, N - 2 ); 
        a = sqrt(larger_eigenval*t0_sqrd);
        b = sqrt(smaller_eigenval*t0_sqrd);
    otherwise
        if strcmp(ellipseType(end-1:end),'CI')
            critVal = str2double(ellipseType(1:end-2))./100;
            if critVal < 1 && critVal > 0                
                % following Eqn. 5-19 of Johnson & Wichern (2007):
                t0_sqrd = ( (N-1)*2 )/( N*(N-2) ) * finv( critVal, 2, N - 2 ); 
                a = sqrt(larger_eigenval*t0_sqrd);
                b = sqrt(smaller_eigenval*t0_sqrd);
            else
                error('CI range must be on the interval (0, 100). Please see the help!')
            end
        else
            error('You entered an invalid error ellipse type. Please see the help!')
        end
end

% the ellipse in x and y coordinates
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
errorEllipse = [ellipse_x_r;ellipse_y_r]' * R;

%Shift to be centered on mean coordinate
errorEllipse = bsxfun(@plus,errorEllipse,meanXy);

% find vector lengths of each point on the ellipse
norms = nan(1,length(errorEllipse));
for pt = 1:length(errorEllipse)
    norms(pt) = norm(errorEllipse(pt,:));
end
[ampMinNorm,ampMinNormIx] = min(norms);
[ampMaxNorm,ampMaxNormIx] = max(norms);
ampEllipseExtremes = [ampMinNorm,ampMaxNorm];
ampDiff = [norm(meanXy) - ampEllipseExtremes(1), ampEllipseExtremes(2) - norm(meanXy)];

% calculate phase angles & find maximum pairwise difference to determine phase bounds
phaseAngles = atan2(errorEllipse(:,2),errorEllipse(:,1));
pairs = combnk(phaseAngles,2); % find all 2-element subsets
diffPhase = abs(pairs(:,2) - pairs(:,1)); % find the absolute difference of each pair
diffPhase(diffPhase > pi) = 2*pi - diffPhase(diffPhase > pi); % unwrap the difference
[~,maxDiffIdx] = max(diffPhase);
anglesOI = pairs(maxDiffIdx,:);
phaseMinIx = find(phaseAngles == anglesOI(1),1);
phaseMaxIx = find(phaseAngles == anglesOI(2),1);

% convert to deg (if necessary) & find difference between (max bound and mean phase) and (mean phase and min bound)
% note: everything converted from [-pi, pi] to [0, 2*pi] for unambiguous calculation
convFactor = [180/pi, 1];
unwrapFactor = [360, 2*pi];
phaseEllipseExtremes = [phaseAngles(phaseMinIx),phaseAngles(phaseMaxIx)]*convFactor(returnRad+1);
phaseEllipseExtremes(phaseEllipseExtremes<0) = phaseEllipseExtremes(phaseEllipseExtremes<0)+unwrapFactor(returnRad+1);
phaseBounds = [min(phaseEllipseExtremes), max(phaseEllipseExtremes)];
meanPhase = atan2(meanXy(2),meanXy(1))*convFactor(returnRad+1);
meanPhase(meanPhase < 0) = meanPhase(meanPhase < 0) + unwrapFactor(returnRad+1);

% if the ellipse overlaps with the origin, defined by whether the phase angles are in all 4 quadrants
phaseAngles(phaseAngles < 0) = phaseAngles(phaseAngles < 0)+2*pi;
if ~isempty(phaseAngles(phaseAngles > 0 & phaseAngles < pi/2)) && ~isempty(phaseAngles(phaseAngles > pi/2 & phaseAngles < pi)) ...
        && ~isempty(phaseAngles(phaseAngles > pi/2 & phaseAngles < 3*pi/2)) && ~isempty(phaseAngles(phaseAngles > 3*pi/2 & phaseAngles < 2*pi))   
    amplBounds = [0, ampMaxNorm];
    maxVals = [360, 2*pi];
    phaseBounds = [0, maxVals(returnRad+1)];
    phaseDiff = [abs(phaseBounds(1) - meanPhase), abs(phaseBounds(2) - meanPhase)];
else
    amplBounds = ampEllipseExtremes;
    phaseDiff = [abs(phaseBounds(1) - meanPhase), abs(phaseBounds(2) - meanPhase)];
end

% unwrap phase diff for any ellipse that overlaps with positive x axis
phaseDiff(phaseDiff > unwrapFactor(returnRad+1)/2) = unwrapFactor(returnRad+1) - phaseDiff(phaseDiff > unwrapFactor(returnRad+1)/2);
    
zSNR = norm(meanXy)/mean([norm(meanXy)-ampMinNorm,ampMaxNorm-norm(meanXy)]);

%% PLOT DATA
if makePlot
    figure;
    hold on;
    plot(xyData(:,srIx),xyData(:,siIx),'ko','MarkerFaceColor','k') 
    axis equal;
    line([0 meanXy(1)],[0 meanXy(2)],'Color','k','LineWidth',1);    
    plot(errorEllipse(:,1), errorEllipse(:,2),'b-','LineWidth',1); 
    hold on;
    plot([meanXy(1) sqrt(smaller_eigenval).*smaller_eigenvec(1)+meanXy(1)],[meanXy(2) sqrt(smaller_eigenval).*smaller_eigenvec(2)+meanXy(2)],'g-','LineWidth',1); 
    plot([meanXy(1) sqrt(larger_eigenval).*larger_eigenvec(1)+meanXy(1)],[meanXy(2) sqrt(larger_eigenval).*larger_eigenvec(2)+meanXy(2)],'m-','LineWidth',1); 
    
    % set plot dimensions
    x_plotmin = min([floor(min(xlim)),-1]);
    x_plotmax = max([ceil(max(xlim)),1]);
    y_plotmin = min([floor(min(ylim)),-1]);
    y_plotmax = max([ceil(max(ylim)),1]);
    xlim([x_plotmin,x_plotmax]);
    ylim([y_plotmin,y_plotmax]);
    
    line([0 0],[min(ylim) max(ylim)],'Color','k')
    line([min(xlim) max(xlim)],[0 0],'Color','k')
    text(.9*min(xlim),.7*min(ylim),[ellipseType ' ellipse'],'FontSize',14,'Color','b');
    text(.9*min(xlim),.6*min(ylim),'larger eigen vec','FontSize',14,'Color','m');
    text(.9*min(xlim),.5*min(ylim),'smaller eigen vec','FontSize',14,'Color','g');
    xlabel('real');
    ylabel('imag');
    drawnow;
    hold off;
    % change figure dimensions
    set(gcf,'units','centimeters');
    fig_pos = get(gcf,'position');
    fig_pos(3) = 15;
    fig_pos(4) = 15;
    set(gcf,'position',fig_pos);
    
    figure;
    hold on;
    axis equal; 
    plot(errorEllipse(:,1), errorEllipse(:,2),'k-','LineWidth',1) 
    plot([0 errorEllipse(ampMinNormIx,1)],[0 errorEllipse(ampMinNormIx,2)],'r:','LineWidth',1);
    plot([0 errorEllipse(ampMaxNormIx,1)],[0 errorEllipse(ampMaxNormIx,2)],'r:','LineWidth',1);
    plot([0 errorEllipse(phaseMinIx,1)],[0 errorEllipse(phaseMinIx,2)],'b:','LineWidth',1);
    plot([0 errorEllipse(phaseMaxIx,1)],[0 errorEllipse(phaseMaxIx,2)],'b:','LineWidth',1);
    
    line([0 meanXy(1)],[0 meanXy(2)],'Color','k','LineWidth',1); % mean vector 
    
    plot([meanXy(1) a.*larger_eigenvec(1)+meanXy(1)],[meanXy(2) a.*larger_eigenvec(2)+meanXy(2)],'m-','LineWidth',1) % half major axis (1)
    plot([meanXy(1) -a.*larger_eigenvec(1)+meanXy(1)],[meanXy(2) -a.*larger_eigenvec(2)+meanXy(2)],'m-','LineWidth',1) % half major axis (2)
    plot([meanXy(1) b.*smaller_eigenvec(1)+meanXy(1)],[meanXy(2) b.*smaller_eigenvec(2)+meanXy(2)],'g-','LineWidth',1) % half minor axis (1)
    plot([meanXy(1) -b.*smaller_eigenvec(1)+meanXy(1)],[meanXy(2) -b.*smaller_eigenvec(2)+meanXy(2)],'g-','LineWidth',1) % half minor axis (2)
    
    text(errorEllipse(ampMinNormIx,1),errorEllipse(ampMinNormIx,2),sprintf('%.2f',ampMinNorm),'FontSize',18,'Color','r')
    text(errorEllipse(ampMaxNormIx,1),errorEllipse(ampMaxNormIx,2),sprintf('%.2f',ampMaxNorm),'FontSize',18,'Color','r')
    text(errorEllipse(phaseMinIx,1),errorEllipse(phaseMinIx,2),sprintf('%.2f',phaseEllipseExtremes(1)),'FontSize',18,'Color','b')
    text(errorEllipse(phaseMaxIx,1),errorEllipse(phaseMaxIx,2),sprintf('%.2f',phaseEllipseExtremes(2)),'FontSize',18,'Color','b')
    
    % set plot dimensions
    x_plotmin = min([floor(min(xlim)),-1]);
    x_plotmax = max([ceil(max(xlim)),1]);
    y_plotmin = min([floor(min(ylim)),-1]);
    y_plotmax = max([ceil(max(ylim)),1]);
    xlim([x_plotmin,x_plotmax]);
    ylim([y_plotmin,y_plotmax]);
    
    text(meanXy(1),meanXy(2),sprintf('%.2f',norm(meanXy)),'FontSize',18)
    text(.9*min(xlim),.7*min(ylim),'ampl. bounds','FontSize',14,'Color','r');
    text(.9*min(xlim),.6*min(ylim),'phase. bounds','FontSize',14,'Color','b');
    text(.9*min(xlim),.5*min(ylim),'mean ampl.','FontSize',14,'Color','k');

    line([0 0],[min(ylim) max(ylim)],'Color','k','LineWidth',1);
    line([min(xlim) max(xlim)],[0 0],'Color','k','LineWidth',1);
    xlabel('real');
    ylabel('imag');
    hold off
    drawnow;
    % change figure dimensions
    set(gcf,'units','centimeters');
    fig_pos = get(gcf,'position');
    fig_pos(3) = 15;
    fig_pos(4) = 15;
    set(gcf,'position',fig_pos);
end

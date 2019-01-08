function varargout=plotewma(ewmaplot,axhan)
%PLOTEWMA  Shows an EWMA control chart for reaction time and accuracy
%   PLOTEWMA(DMATEWMAPLOT), shows a control chart. DMATEWMAPLOT is a
%   structure obtained from EWMAV2.
%
%   PLOTEWMA(DMATEWMAPLOT, AXHAN), places the control chart in the axes
%   with handle AXHAN.
%
%   H = PLOTEWMA(...), returns a handle to the axes where the control
%   chart was placed.
%
%   See also EWMAV2, OUTLIERTREATMENT.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.3: Added option for switching axes.

%% Input checking
if nargin<1
    error('DMAT:plotewma:notEnoughInputs',...
        'EWMAPLOT requires at least one input argument.')
end
if nargin==2
    prevax = gca; % Store current axes
    axes(axhan);  % Make selected axes active
end

if isempty(ewmaplot)
    disp('Input structure was empty!')
    return
end

%% Create a patch inside the control limits
h=fill([ewmaplot.RTs ewmaplot.RTs(end:-1:1)],...
    [ewmaplot.UCL ewmaplot.LCL],'r');
set(h,'FaceAlpha',.1,'EdgeColor',[.95 .55 .55])

%% Make the plot area prettier
ax=gca;
set(ax,'XLim',[0 1.4])
ylabel('State of the system (c)')
xlabel('Reaction time')
title(sprintf('EWMA control chart (\\lambda=%g, L=%g)',...
    ewmaplot.Lambda,ewmaplot.Limits))

%% Plot the EWMA statistic in red and blue
% If less than 500 data points, plot dots, otherwise, lines.
if length(ewmaplot.BlueX)+length(ewmaplot.RedX)<500
    bf = 'b.'; rf = 'r.';
else
    bf = 'b-'; rf = 'r-';
end
hold on
plot(ewmaplot.BlueX,ewmaplot.BlueY,bf);
plot(ewmaplot.RedX,ewmaplot.RedY,rf);

%% Also plot the control state
line([ewmaplot.RTs(1) 1.4],ewmaplot.ControlState([1 1]),...
    'Color',[0 0 0],'LineStyle',':')
hold off

%% If requested, return axis handle
if nargout
    varargout{1}=ax;
end
if nargin==2 % Put GCA back to what it was
    axes(prevax);
end
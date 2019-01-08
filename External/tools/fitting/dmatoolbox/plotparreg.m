function varargout=plotparreg(output,varargin)
%PLOTPARREG  Plot deviance function around one or two given parameters
%   PLOTPARREG(OUTPUT,P), where P is a scalar indexing the free parameters
%   in the model used in OUTPUT, shows a 2D plot of the deviance function
%   around the minimum.
%
%   PLOTPARREG(OUTPUT,P1,P2), where P1 and P2 are both indexing scalars,
%   shows a 3D (surface) plot.
%
%   F = PLOTPARREG(...) returns a vector or matrix of the parameter space
%   to plot.
%
%   The scalar indexing arguments can also be strings containing an index,
%   or can be strings describing the parameters, like 'a(1)' for the first
%   element of the design matrix for a.
%
%   See also FITLAST, NAMEPARS.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%% Check input
if nargin<2
    error('DMAT:plotparreg:notEnoughInputs',...
        'PLOTPARREG requires at least two input arguments.')
end
if numel(output)~=1
    error('DMAT:plotparreg:multipleOutputStructures',...
        'Argument OUTPUT should be a single structure.')
end


%% Basic settings
N1 = 25;          % columns
N2 = 25;          % rows
Nsingle = 50;
offset1 = 1e-2;
offset2 = 1e-2;


%% Identify used parameters
allp = namepars(output);
npar = length(output.Options.controls.small);

if any([varargin{2:end}]>npar) || any([varargin{2:end}]<1)
    error('DMAT:plotparreg:parameterIndexTooHigh',...
        'At least one parameter index doesn''t exist.')
end

%% Process input
p = nan(1,nargin-1);
for c1=1:nargin-1
    for c2=1:npar
        up = [allp{c2,1} '(' num2str(allp{c2,2}) ')'];
        if strcmp(varargin{c1},up)
            p(c1)=c2;
        elseif all(isstrprop(varargin{c1},'digit'))
            p(c1)=str2double(varargin{c1});
        elseif isnumeric(varargin{c1})
            p(c1)=varargin{c1};
        end
    end
end
if any(isnan(p))
    error('DMAT:plotparreg:BadInput','Bad input to PLOTPARREG')
end
pna = allp(p,:);

%% Name parameters
for a=1:length(p)
    fprintf('   Parameter %u is %s(%i):%8.5f.\n',p(a),pna{a,1},pna{a,2},pna{a,3});
end

%% Online plotting?
online = nargout==0;

%% Plotting
% one dimension or two?
vin = 2;
if length(p)<2,p(2)=p;end
if p(1)==p(2),N1=1;N2=Nsingle;vin=[90 0];grid off;end

%% Define domain
min1 = output.Options.controls.small(p(1));
min2 = output.Options.controls.small(p(2));
range1 = sort([linspace(min1-offset1,min1+offset1,N1-1) min1]);
range2 = sort([linspace(min2-offset2,min2+offset2,N2-1) min2]);
[r1 r2] = meshgrid(range1,range2);

%% Initialise other parameters
par = output.Options.controls.small;
of = output.Options.objecfun;
X2p = NaN(N2,N1);

%% Pinpoint recovered minimum
if online,
    plot3(min1,min2,0,'ro','MarkerSize',5)
    hold on
end

%% Calculate all points
for c1=1:N1
    for c2=1:N2
        par(p) = [range1(c1) range2(c2)];
        cc = of(par)-output.Fitvalue;
        if cc>1e3, cc=NaN; end
        X2p(c2,c1) = cc;
        if N1==1,X2p(c2,2)=cc;end
        if online
            surf(r1,r2,X2p,'EdgeColor','interp')
            view(vin)
            drawnow
        end
    end
end

%% Put minimum again
if online
    hold on
    plot3(min1,min2,0,'ro','MarkerSize',5)
    hold off
    axis tight
    try
        xlabel(pna(1))
        ylabel(pna(2))
    catch
        title(pna)
    end
end

%% Found any points below minimum?
if any(X2p<-output.Options.ObjectiveDecimals)
    fprintf('   Found a better point: X? = %17.16f.\n',output.Fitvalue+min(X2p(:)));
    fprintf('   Compare to recovered: X? = %17.16f.\n',output.Fitvalue);
end

%% Output
if ~online
    varargout{1}=X2p;
end
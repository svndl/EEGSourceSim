function varargout = edfcdf(data,N,rep,param,quantiles,handle)
%EDFCDF Plot cumulative distributions of reaction time
%   EDFCDF(RTS,CONDN,REP,PARAM,BINEDGES,HANDLE), plots an empirical and
%   theoretical cumulative distribution of reaction time, for the reaction
%   times in the vector RTS, for reply REP, overlaying a theoretical
%   cumulative distribution with parameters PARAM, adding reference lines
%   located at the values in BINEDGES, in the axes indicated by HANDLE.
%   RTS should contain only the reaction times for the specific reply,
%   while N should contain the number of data points in this condition
%   i.e., both replies).
%
%   EDFCDF(RTS,CONDN,REP,PARAM,BINEDGES), makes a new figure.
%
%   H = EDFCDF(...), returns a handle to the axes.
%
%   If the fourth and/or fifth argument are the empty matrix, theoretical
%   distributions and/or reference lines are omitted (but note that you
%   have to provide both arguments regardless).
%
%   See also QPPLOT.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Made plot prettier.
%  Edit 0.4: Renamed variable tcdf to theocdf because MATLAB was being
%  silly about it.

%% Check input
if nargin<6
    if nargin<5
        error('DMAT:edfcdf:NotEnoughInputs',...
            'EDFCDF requires at least five input arguments.')
    end
    figure
    handle = axes;
end

set(handle,'NextPlot','add','Box','on')

if rep
    yl = 'P(X=1)'; % label
    col = [0 1 0]; % plot color    
else
    yl = 'P(X=0)';
    col = [1 0 0];
end


%% Empirical cumulative distribution
cn = size(data,1);
maxrt = 3;
data(data>maxrt)=maxrt;
if cn>0
    edom = sort(data);
    maxrt = max(data);
    ecdf = (1:cn)./N;
else
    edom = [0 maxrt];
    ecdf = [0 0];
end

%% EDF
axes(handle)
plot(handle,edom,ecdf,'Color',[0 0 1])

%% Theoretical CDF
if ~isempty(param)
    tdom = linspace(0,maxrt,30);
    theocdf = cdfdif(tdom,rep,param);
    plot(handle,tdom,theocdf,'Color',col)
end

%% Lines
if ~isempty(quantiles)
    ns = length(quantiles);
    line(quantiles([1 1],:),[zeros(ns,1) ones(ns,1)]',...
        'LineStyle',':','Color',[.55 .55 .55])
end

%% Pretty
xlabel(handle,'Time (seconds)','FontSize',8)
ylabel(handle,yl,'FontSize',8)
axis([0 maxrt 0 max(theocdf(end),cn/N)])

%% Finalize output
if nargout
    varargout{1} = handle;
end
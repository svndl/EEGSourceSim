function dcel = processdata(data,perc,edges)
%PROCESSDATA  Calculates statistics on each condition
%   DCEL = PROCESSDATA(DATA,PERC,EDGES), where DATA is an N-by-2 data
%   matrix (for one condition) of the format [response seconds], returns
%   DCEL, an 8-by-1 cell matrix with summary statistics for that condition,
%   as described in the documentation for FITDIFFV13. PERC may contain
%   percentiles for determining bin edges (data-dependent bin edges), or
%   may be empty. EDGES may contain bin edges (a-priori edges), or may be
%   empty. Both vectors may be omitted or may be empty, in which case
%   default values are used. Default values are also used if both are
%   nonempty.
%
%   See also FITDIFFV13, SPLITDATA.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%% Check input
if nargin<3
    edges = [];
    if nargin<2
        perc = [];
        if nargin<1
            error('DMAT:processdata:notEnoughInputs',...
                'PROCESSDATA requires at least one input.')
        end
    end
end

%% Check if DATA not a cell matrix of split data
if iscell(data)
    dcel = cell(size(data));
    for ctr = 1:length(data)
        dcel(ctr) = {processdata(data{ctr},perc,edges)};
    end
    return
end

%% Split dataset
yeses = logical(data(:,1));
nos = sort(data(~yeses,2));
yess = sort(data(yeses,2));
nyes = sum(yeses);
nno = sum(~yeses);
total_n = nno+nyes;
nO = [];
yO = [];
nQ = [];
yQ = [];

if ~isempty(edges) && isempty(perc)
%% User-defined bins
    %disp 'FixedBinEdges: User-defined'
    yQ = edges(1,:);
    nQ = edges(end,:);
elseif isempty(edges) && ~isempty(perc)
%% Find bin edges from user-defined percentiles
    %disp 'FixedBinEdges: Find from user-defined percentiles'
    yP = perc(1,:);
    if nyes>=11;
        spc = [0 100*(0.5:(nyes-0.5))./nyes 100]';
        vec = [yess(1,:); yess(1:nyes,:); yess(nyes,:)];
        yQ = interp1q(spc,vec,yP')';
    end
    nP = perc(end,:);
    if nno>=11;
        spc = [0 100*(0.5:(nno-0.5))./nno 100]';
        vec = [nos(1,:); nos(1:nno,:); nos(nno,:)];
        nQ = interp1q(spc,vec,nP')';
    end
else
%% Default bin edges
    %disp 'FixedBinEdges: Default'
    yQ = [0.30 0.36 0.42 0.52 0.80];
    nQ = [0.38 0.47 0.56 0.70 1.00];
end

if nyes>=11;
    yO = histc(yess,[0 yQ inf])';
    yO = yO(1:end-1);
end
if nno>=11;
    nO = histc(nos,[0 nQ inf])';
    nO = nO(1:end-1);
end

%% Find extreme reaction times for adaptive filter
extrt = [min([nos;yess]),max([nos;yess])];
dcel = {yQ;nQ;yO;nO;nyes;nno;total_n;extrt};
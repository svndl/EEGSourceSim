function varargout=qpplot(data,perc,param,handle)
%QPPLOT  Quantile probability plots
%   QPPLOT(DATA,PERC,PARAM,HANDLE), displays a quantile probability plot,
%   for the N-by-3 data in DATA, using bins based on percentiles in PERC,
%   overlaying theoretical values from parameter set PARAM, in the axes
%   indicated by HANDLE. Note that PERC are percentiles from 0 to 100.
%
%   QPPLOT(DATA,PERC,PARAM), makes a new figure.
%
%   QPPLOT(DATA,PERC,[]) or QPP(DATA,PERC,[],HANDLE) do not plot the model
%   fit lines (but note that you have to provide a third argument
%   regardless).
%
%   H = QPPLOT(...), returns a handle to the axes.
%
%   QPPLOT(OUTPUT) is a quick-access version that only plots model
%   predictions in a new figure, with PERC = [10 30 50 70 90]; This version
%   assumes that conditions are sorted by difficulty, with the easiest
%   coming first.
%
%   See also EDFCDF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.4: Edited use of INVCDFDIF, which no longer takes defective
%   percentiles as input.
%   Edit 0.4: Moved PCTILE to own M-file.
%   Edit 0.4: Added quick-access version.

%% Check input
if nargin<4
    if nargin<3
        if nargin==1 && isstruct(data) && numel(data)==1 % Quick plot!
            if ~isfield(data,'FitInfo')||~isfield(data.FitInfo,'PredictedQuantiles')
                error('DMAT:qpplot:predictionsFieldNotFound',...
                    'The FitInfo structure does not have a PredictedQuantiles field (this is new to DMAT 0.4).')
            end
            plot(1:data.Options.nsets*2,...
                data.FitInfo.PredictedQuantiles([1:2:data.Options.nsets*2 2*data.Options.nsets:-2:2],:)',...
                '.-')
            ylabel('Time (seconds)','FontSize',8)
            xlabel('Condition (split by response and performance)','FontSize',8)
            grid on
            return
        else
        error('DMAT:qpplot:NotEnoughInputs',...
            'QPPLOT requires three input arguments, or one output structure.')
        end
    end
    figure
    handle = axes;
end
if ~isvector(perc)
    error('DMAT:qpplot:percentileVectorNotVector',...
        'Input variable PERC should be a vector.')
end

if all(perc<=1)
    warning('DMAT:qpplot:SmallPercentiles',...
        ['QPPLOT requires that the PERC input is in the range 0-100, ' ...
        'but yours are in the range 0-1.'])
end
if any(perc<=0|perc>=100)
    error('DMAT:qpplot:InvalidPercentiles',...
        ['QPPLOT requires that the PERC input is in the range 0-100 ' ...
        '(exclusive).'])
end
if ~ishandle(handle)
    error('DMAT:qpplot:InvalidHandle',...
        'Input HANDLE does not seem to be a valid handle.')
elseif ~strcmp(get(handle,'Type'),'axes')
    error('DMAT:qpplot:NotAnAxisHandle',...
        'Input HANDLE does not seem to be a handle to an AXES object.')
end
if ~isvaliddataset(data)
    error('DMAT:qpplot:InvalidDataset',...
        'Input DATA does not seem to be a valid dataset.')
end

%% Determine accuracies
conds = unique(data(:,1));
nsets = length(conds);
for ctr=1:nsets
    this_cond = data(data(:,1)==conds(ctr),2:3);
    yeses = (this_cond(:,1)==1);
    percentage_correct(1,ctr) = mean(yeses);
    rt_yeses(:,ctr) = pctile(this_cond(yeses,2),perc)';
    rt_noes(:,ctr) = pctile(this_cond(~yeses,2),perc)';
end
percentage_incorrect = 1-percentage_correct;

%% Build plot
abscissa = [percentage_incorrect percentage_correct(end:-1:1)];
ordinate = [rt_noes rt_yeses(:,end:-1:1)];
[abscissa index]=unique(abscissa);
ordinate = ordinate(:,index);

%% Plot it
plot(handle,abscissa,ordinate','x')

%% Plot model recovery if able
rows = size(param,1);
ordin2 = [];
if ~isempty(param)
    if rows~=nsets
        error('DMAT:qpplot:DataAndParametersInconsistent',...
            ['The dataset in DATA contains %i conditions, but ',...
            'there are (is) %i parameter vector(s) in PARAM'],nsets,rows)
    end
    ordin2 = zeros(length(perc),rows*2);
    for ctr=1:rows
        ordin2(:,ctr) = invcdfdif(perc/100,0,param(ctr,:))';
        ordin2(:,end+1-ctr) = invcdfdif(perc/100,1,param(ctr,:))';
    end
    ordin2 = ordin2(:,index);
    hold on
    plot(abscissa,ordin2','-');
    hold off
end

ll = length(index);
index(index>ll/2) = ll+1-index(index>ll/2);
xticklab = [cellfun(@(x)sprintf('W%i',x),num2cell(index(1:end/2)),'uni',0),...
    cellfun(@(x)sprintf('C%i',x),num2cell(index(end/2+1:end)),'uni',0)];

%% Finalize plot
set(handle,...
    'XTick',abscissa,...
    'XTickLabel',xticklab,...    'XLim',[0 1],...
    'YLim',[min([ordinate(:);ordin2(:)])*.8 max([ordinate(:);ordin2(:)])*1.1]);
ylabel('Time (seconds)','FontSize',8)
xlabel('Condition (split by response and sorted by performance)','FontSize',8)

%% And return output if requested
if nargout
    varargout{1}=handle;
end
function PlotScalp(pattern,title_text)
    % simplifed function to just plot scalp pattern
    if ~exist('title_text','var')
        title_text = [] ;
    end
    
    load(fullfile('Electrodeposition.mat')); tEpos =tEpos.xy;% electrode positions used for plots

    colorbarLimits = [min(pattern),max(pattern)];
    conMap = jmaColors('coolhotcortex');
    Probs{1} = {'facecolor','none','edgecolor','none','markersize',10,'marker','o','markerfacecolor','g' ,'MarkerEdgeColor','k','LineWidth',.5};% plotting parameters

    mrC.plotOnEgi(pattern,colorbarLimits,false,[],false,Probs); 

    if ~isempty(title_text)
        title(title_text)
    end
function printpdf(h,savePath,format)
    % first use the same non-relative unit system for paper and screen (see
    % below)
    set(h,'PaperUnits','centimeters');

    % now get all existing plots/subplots
    a=get(h,'Children');
    nfigs=length(a);

    % bounds will contain lower-left and upper-right corners of plots plus one
    % line to make sure single plots work
    bounds=zeros(nfigs+1,4);
    bounds(end,1:2)=inf;
    bounds(end,3:4)=-inf;

    % generate all coordinates of corners of graphs and store them in
    % bounds as [lower-left-x lower-left-y upper-right-x upper-right-y] in
    % the same unit system as paper (centimeters here)
    for i=1:nfigs
        set(a(i),'Unit','centimeters');
        pos=get(a(i),'Position');
        inset=get(a(i),'TightInset');
        bounds(i,:)=[pos(1)-inset(1) pos(2)-inset(2) ...
            pos(1)+pos(3)+inset(3) pos(2)+pos(4)+inset(4)];
    end

    % compute the rectangular convex hull of all plots and store that info
    % in mypos as [lower-left-x lower-left-y width height] in centimeters
    auxmin=min(bounds(:,1:2))-2;
    auxmax=max(bounds(:,3:4))+2;
    mypos=[auxmin auxmax-auxmin];

    % set the paper to the exact size of the on-screen figure using
    % figure property PaperSize [width height]
    set(h,'PaperSize',[mypos(3) mypos(4)]);

    % ensure that paper position mode is in manual in order for the
    % printer driver to honor the figure properties
    set(h,'PaperPositionMode', 'manual');

    % use the PaperPosition four-element vector [left, bottom, width, height]
    % to control the location on printed page; place it using horizontal and
    % vertical negative offsets equal to the lower-left coordinates of the
    % rectangular convex hull of the plot, and increase the size of the figure
    % accordingly
    set(h,'PaperPosition',[-mypos(1) -mypos(2) ...
        mypos(3)+mypos(1) mypos(4)+mypos(2)]);

    % print stuff
    if strcmp(format,'png')
        print(gcf,'-r600','-dpng',[savePath,'png']) 
    else
        print('-dpdf',[savePath,'pdf']);
    end
end
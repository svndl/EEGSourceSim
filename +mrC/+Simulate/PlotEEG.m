function h = PlotEEG(ASDEEG,Freq,savepath,subID, masterList,signalFF,SignalType,jcolors,Mode,EOI,FOI,A,W)
% This function provides a dual plot of electrode amplitude/phase spectrum and topographic
% map of amplitude/phase at a specific frequency (similar to powerDiva)

% This function is interactive: click on the head plot to select electrode
% and click on spectrum plot to select frequency bin

% PRESS ENTER: to save the current figure in the save folder path
% PRESS ESC: to exit the plot
% PRESS N to change the normalization of the head plot

% INPUT:
    % ASDEEG: nf x nelec matrix, amplitude spectrum of EEG, nf is number of
            % frequency bins and nelec is number of electrodes
    % Freq: nf x 1 vector, indicating the frequency bins
    % savepath: the string input indicating the folder to save the plots
    % master list: a cell array containing the names of ROIs
    % signalFF: fundamental frequency of the seed signals
    % SignalType: plot amplitude or phase, [Amplitude]/Phase
    % Mode: the mode of plotting, can be interactive or simple, in
            % interactive mode user can select electrodes and frequency bin
            % [Interact]/Simple
    % EOI: the electrode to plot the (initial) spectrum
    % FOI: the frequency to plot the (initial) spectrum topo map
    
    
% Written by Elham Barzegaran,3.7.2018
% modified by EB: 6.6.2018
% modified by sebastian bosse 8.8.2018
% modified by EB: 8.21.2018
% modified by sb: 9.17.2018

%% set parameters
if ~exist('SignalType','var')|| isempty(SignalType) % plot phase or amp
    if min(ASDEEG(:))>0
        SignalType = 'Amplitude';
    else
        SignalType = 'Phase';
    end
end

FS=14;%font size
Fmax = numel(Freq); % maximum frequency

if ~exist('EOI','var') || isempty(EOI)
    EOI = 83; % electrode used in the plot
end

load('Electrodeposition.mat'); tEpos =tEpos.xy;% electrode positions used for plots

if ~exist('FOI','var') || isempty(FOI)
    if isempty(signalFF), 
        FFI=2; 
    else 
        FFI = signalFF;
    end
else
    FFI = FOI;
end
[~,~,FFI] = intersect(FFI,round(Freq*1000)/1000);% find the index of fundamental frequencies in frequency bins and keep them for plotting purpose
FOI = FFI(1);

if ~exist('Mode','var') || isempty(Mode)
   Mode = 'Interact';
end

% allow for visualization of spatial filters
if ~exist('A','var') || isempty(Mode)
    space = 'chan';
else
    space = 'comp';
    EOI = 1 ; % start with the first component
    A = abs(A) ; % visualization in amplitude domain
    %Mode = '';
    weight = 1;
end

%% Plot prepration
Probs{1} = {'facecolor','none','edgecolor','none','markersize',10,'marker','o','markerfacecolor','g' ,'MarkerEdgeColor','k','LineWidth',.5};% plotting parameters
if ~exist('jcolors','var') || isempty(jcolors)
    switch SignalType
        case 'Amplitude'
            conMap = jmaColors('hotcortex');
        case 'Phase'
            conMap = jmaColors('phasecolor');
            %conMap = hsv(64);
    end
else
    conMap = jmaColors(jcolors);
end

h=figure;
set(h,'units','centimeters')
set(h, 'Position',[1 1 35 16]);
set(h,'PaperPositionMode','manual')


subplot(2,2,2),axis off;% Show simulated signal information

nrs = max(numel(masterList),3);
text(.1,1,['' subID],'fontsize',FS+1,'fontweight','bold');
for i= 1:numel(masterList)
    if ~isempty(signalFF)
        text (.1,1-((.15)*(i)),['Source' num2str(i) ': f' num2str(i) ' = ' num2str(signalFF(i)) ' Hz,  ' strrep(masterList{i},'_','-')],'fontsize',FS-2);
    else
        text (.1,1-((.15)*(i)),['Source' num2str(i) ': ' strrep(masterList{i},'_','-')],'fontsize',FS-2);
    end
end
set(gca,'tag','info');


%% draw buttons for component selection
if strcmpi(space,'comp')
    ax1 = axes('Units','pixel','Position',[30 400 50 10],'Visible','off');
    axes(ax1);
    text(0, 0,'RC number','fontweight','bold')
    popup = uicontrol('Style', 'popup',...
                      'String', '1|2|3|4|5|6|7|8|9',...
                      'Position', [20 340 100 50],...
                      'Callback',@set_component);
    % plot which weight
    ax2 = axes('Units','pixel','Position',[30 360 50 10],'Visible','off');
    axes(ax2);
    text(0, 0,'Weight','fontweight','bold')
    popup = uicontrol('Style', 'popup',...
                  'String', 'A|W',...
                  'Position', [20 300 100 50],...
                  'Callback',@set_weight);
end
N = 1;
colorbarLimits = [] ;



setColorbar()
%---------------------------HEAD PLOT----------------------------------
plot_topog();
%---------------------------SPECTRUM PLOT------------------------------
plot_spectrum();
%----------------------------------------------------------------------

%% ----------------Reads keyboard or mouse click-----------------------
if strcmpi(Mode,'Interact')
    set(h, 'WindowKeyPressFcn',@keyPressCallback);
    set(h, 'WindowButtonDownFcn', @buttonPressCallback);
end


%----------------------------------------------------------------------
    function set_component(source,event)
        EOI = get(source,'value');
        plot_topog();
        plot_spectrum();
    end

    function set_weight(source,event)
        weight = get(source,'value');
        plot_topog();
    end
    
    function keyPressCallback(source, eventdata)
        key = eventdata.Key;
        if key==27 % (Esc key) -> close the plot and return
            close;
            h=[];
            return
        elseif key==13 % (Enter key) -> save the current figure
            set(h, 'PaperPosition',[1 1 12 5]);
            print(fullfile(savepath,['SimEEG_Subject' subID 'Electrode' num2str(EOI) '_Freq' num2str(Freq(FOI)) 'Hz.tif']),'-dtiff','-r300');% Later I can update this to contain the simulation parameters
        elseif strcmp(key,'n')||strcmp(key,'N') % if n or N is pressed -> change head plot normalization
            if N==1, N=0; else N=1;end

        elseif strcmpi(key,'leftarrow') && strcmpi(space,'comp')  % cursor left
            EOI = mod(EOI-1,size(A,2)+1) ;
            if EOI==0 % 0 is not valid
                 EOI = mod(EOI-1,size(A,2)+1) ;
            end
        elseif strcmpi(key,'rightarrow') && strcmpi(space,'comp') % cursor right
            EOI = mod(EOI+1,size(A,2)+1) ;
            if EOI==0 % 0 is not valid
            	EOI = mod(EOI+1,size(A,2)+1) ;
            end
        elseif strcmpi(key,'uparrow') && strcmpi(space,'comp')  % cursor up
            EOI = mod(EOI+10,size(A,2)+1) ;
            if EOI==0 % 0 is not valid
              EOI = mod(EOI+1,size(A,2)+1) ;
            end
        elseif strcmpi(key,'downarrow') && strcmpi(space,'comp')  % cursor down
            EOI = mod(EOI-10,size(A,2)+1) ;
            if EOI==0 % 0 is not valid
                EOI = mod(EOI-1,size(A,2)+1) ;
            end
        elseif key=='a' 
            FOI = mod(FOI-1,Fmax+1) ;
            if FOI==0 % 0 is not valid
                FOI = mod(FOI-1,Fmax+1) ;
            end
        elseif key=='d' 
            FOI = mod(FOI+1,Fmax+1) ;
            if FOI==0 % 0 is not valid
                FOI = mod(FOI+1,Fmax+1) ;
            end
        elseif key=='w' 
            FOI = mod(FOI+10,Fmax+1) ;
            if FOI==0 % 0 is not valid
                FOI = mod(FOI+1,Fmax+1) ;
            end
        elseif key=='s' 
          FOI = mod(FOI-10,Fmax+1) ;
          if FOI==0 % 0 is not valid
                FOI = mod(FOI-1,Fmax+1) ;
          end
        end
        % update plots
        setColorbar();
        plot_topog();
        plot_spectrum();
    end
    
    function buttonPressCallback(source, eventdata)
        mousept = get(gca,'currentPoint');
        SPI = get(gca,'tag');
        x = mousept(1,1);
        y = mousept(1,2);
        % update location
        switch SPI

        case '1'
            if ~strcmpi(space,'comp')
                Epos2= repmat([x y],[128 1]);
                dis = sqrt(sum((tEpos-Epos2).^2,2));
                [~,EOI] = min(dis);
            end
        case '2'
            [~,FOI] = min(abs(repmat(x,[1 size(ASDEEG,1)])-Freq));
        end
        % update plots
        setColorbar();
        plot_topog();
        plot_spectrum();
    end
    
%--------------------------------------------------------------------------
   function setColorbar()
    if strcmpi(SignalType,'Amplitude')
        if strcmpi(space,'comp') 
            colorbarLimits = [-0 max(A(:,EOI)*ASDEEG(FOI,EOI)')];
        else
            if N == 1, colorbarLimits = [-0 max(ASDEEG(FOI,:))];
            else colorbarLimits = [-0 max(ASDEEG(:))];
            end
        end
    elseif strcmpi(SignalType,'Phase')
        colorbarLimits = [min(ASDEEG(:)) max(ASDEEG(:))];
    end
   end

%--------------------------------------------------------------------------
    function plot_topog()
        % topography map plot
        
        sp1 = subplot(1,2,1);delete(sp1);
        sp1 = subplot(1,2,1);
        if strcmpi(Mode,'Interact') && ~strcmpi(space,'comp'),
            mrC.plotOnEgi(ASDEEG(FOI,:)',colorbarLimits,false,EOI,false,Probs); 
        else
            if strcmpi(space,'comp') 
                if weight==1
                    conMap = jmaColors('hotcortex');
                    mrC.plotOnEgi(A(:,EOI)*ASDEEG(FOI,EOI)',colorbarLimits,false,[],false,Probs); 

                else
                    mrC.plotOnEgi(W(:,EOI),[-max(abs(W(:,EOI))) max(abs(W(:,EOI)))],false,[],false,Probs);
                    conMap = jmaColors('coolhotcortex');
                end
            else
                mrC.plotOnEgi(ASDEEG(FOI,:)',colorbarLimits,false,[],false,Probs); 
            end
        end
        title(['Frequency = ' num2str(Freq(FOI)) 'Hz'],'fontsize',FS);
        set(sp1,'tag',num2str(1));
        colormap(conMap);
        colorbar;
    end
%--------------------------------------------------------------------------
    function plot_spectrum() 
         % spectrum plot
        SPEC = ASDEEG(1:Fmax,EOI);
        sp2 = subplot(2,2,4); delete(sp2);
        sp2 = subplot(2,2,4); 
        bar(Freq(1:Fmax), SPEC,.15); 
        xlim([Freq(1) Freq(Fmax)]);
        xlabel('Frequency(Hz)','fontsize',FS-2);
        if strcmp(SignalType,'Amplitude')
            ylim([min(SPEC) max(SPEC)*1.1]);
            ylabel('ASD','fontsize',FS-2);
        elseif strcmp(SignalType,'Phase')
            ylim([min(SPEC)*1.1 max(SPEC)*1.1]);
            ylabel('Phase(degree)','fontsize',FS-2);
        end
        set(sp2,'tag',num2str(2));
        
        if strcmpi(Mode ,'Interact'),hold on;  bar(Freq(FOI), SPEC(FOI),.4,'FaceColor','g','EdgeColor','g'); end

        if strcmpi(space,'comp')
            title(['Component ' num2str(num2str(EOI))],'fontsize',FS);
        else
            title(['Electrode ' num2str(num2str(EOI))],'fontsize',FS);
        end
    end
end
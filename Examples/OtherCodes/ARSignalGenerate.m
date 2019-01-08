function ARSignalGenerate(dorun,dosave)
    FigPath = 'Figures';
    addpath(genpath('../BrainNetSimulation'));
    
    if ~exist('dorun','var') || isempty(dorun)
        dorun = false;
    end
    
    if ~exist('dosave','var') || isempty(dosave)
        dosave = false;
    end
    %% design connected network
    if dorun
        Net1 = BrainNetSim(3,300,.98);
        % define net dynamics
        NS = 100000; % length of simulation: time points
        Net1 = Net1.AddNodeFreqs([1],{[8 20]});
        Net1 = Net1.AddNode('AR');
        Net1 = Net1.AddNodeFreqs([4],{[8 20]});
        
        % Add connections to the network
        Net1 = Net1.AddConnection([1 2],'Type','low','LF',13,'Order',24,'Gain',1);
        Net1 = Net1.AddConnection([2 3],'Type','delay','Order',15,'Gain',2);
        Net1 = Net1.AddConnection([1 3],'Type','high','HF',13,'Order',24,'Gain',1/4);

        Net1 = Net1.GenerateARMatrix;
        [Net1,TS_connect] = Net1.Realization(NS);
        Name = 'SimSignals_connect_3Nodes-v2';
        plot_net_sig(Net1,TS_connect,0,[],Name);
        Net_connect = Net1;

        %% design unconnected network
        Net2 = BrainNetSim(3,300,.98);
        % define net dynamics
        NS = 100000; % length of simulation: time points
        Net2 = Net2.AddNodeFreqs([1 2 3],{[8 20],[8],[8 20]});
        Net2 = Net2.AddNode('AR');
        Net2 = Net2.AddNodeFreqs([4],{[8 20]});
        
        Net2 = Net2.GenerateARMatrix;
        [Net2,TS_unconnect] = Net2.Realization(NS);
        Name = 'SimSignals_unconnect_3Nodes-v2';
        plot_net_sig(Net2,TS_unconnect,0,[],Name);
        Net_unconnect = Net2;
        if dosave
            save(fullfile(pwd,'NetTimeSeries_3Nodes-v2'),'Net_unconnect','TS_unconnect','Net_connect','TS_connect');
        end
    else
        load(fullfile(pwd,'NetTimeSeries_3Nodes'));
        
        Name = 'SimSignals_connect_3Nodes';
        plot_net_sig(Net_connect,TS_connect,1,[],Name);
        
        Name = 'SimSignals_unconnect_3Nodes';
        plot_net_sig(Net_unconnect,TS_unconnect,1,[],Name);
    end
end


%%

function plot_net_sig(Net,TS,issave, FigPath,Name)
    % Plots Network realizations
    if ~exist('issave','var') || isempty(issave)
        issave = 0;
    end
    if ~exist('FigPath','var') || isempty(FigPath)
        if ~exist('Figures','dir')
            mkdir('Figures');
        end
        FigPath = fullfile(pwd,'Figures');
    end
    if ~exist('Name','var') || isempty(Name)
        Name = 'NetResults';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fig1 = figure;
    load(fullfile('ExampleData','ROI_colors_Poster.mat'));
    FS = 22;
    for node = 1:Net.NodeNum
        % 
        TSP = TS(node,1:Net.SF*2);
        subplot(Net.NodeNum,2,(node-1)*2+1),plot(TSP,'color','k');
        set(gca,'xtick',0:Net.SF:Net.SF*2,'xticklabel',0:2, 'yticklabel',[],'fontsize',FS);
        xlim([0 Net.SF*2])
        if node == Net.NodeNum, xlabel('time (s)','fontsize',FS); end
        ylabel(['S' num2str(node)],'fontweight','bold','fontsize',FS,'color',Colors(node,:));
        if node ==1 
            T1 =  title('Temporal dynamic','fontsize',FS);
            %set(T1,'position',get(T1,'position')+[0 0.1 0]);
        end
        ylim([-max(abs(TSP)) max(abs(TSP))]*1.1);

        subplot(Net.NodeNum,2,(node-1)*2+2);
        [Z,f] = pwelch(TS(node,:),Net.SF,[],[],Net.SF);
        plot(f,Z,'linewidth',2,'color','k');
        set(gca, 'yticklabel',[],'fontsize',FS-2);
        if node==1
            T2 = title('Power Spectrum Density','fontsize',FS);
            %set(T2,'position',get(T2,'position')+[0 0.1 0]);
        end
        if node == Net.NodeNum,  xlabel('Frequency(Hz)','fontsize',FS);end
        xlim([0 30]);
        ylim([0 max(Z)*1.1]);
        %set(gca,'yticklabel',[])
        
    end
    if issave
        set(Fig1,'PaperPosition',[1 1 8 6])
        print(fullfile(FigPath,Name),'-r300','-dtiff')
    end
end


% Generate ARAX signal
clear; clc;
close all

FigPath = 'C:\Users\Elhamkhanom\Documents\My works\StanfordWorks\Simulation\papers\First paper\Figures\';

%% design network
Net1 = BrainNetSim(3,100,.96);
Net1 = Net1.AddNode('AR');
% define net dynamics
NS = 20000; % length of simulation: time points
Net1 = Net1.AddNodeFreqs([1 2 3],{[9],[4 15],[5 12]});

Net1 = Net1.GenerateARMatrix;
[Net1,TS] = Net1.Realization(NS);

Name = 'SimSignals_Unconnect';
plot_net_sig(Net1,TS,1,FigPath,Name);

%% Add connections to the network

Net1 = Net1.AddConnection([1 2],'Type','bandpass','LF',5,'HF',15,'Order',6,'Gain',1/5);
Net1 = Net1.AddConnection([3 4],'Type','high','HF',7,'Gain',1/4);
Net1 = Net1.AddConnection([1 4],'Type','delay','Order',5,'Gain',1/3);

Net1 = Net1.GenerateARMatrix;
[Net1,TS] = Net1.Realization(NS);

Name = 'SimSignals_connect';
plot_net_sig(Net1,TS,1,FigPath,Name);

%%
function plot_net_sig(Net,TS,issave, FigPath,Name)
    % Plots Network realization
    if ~exist('issave','var') || isempty(issave)
        issave = 0;
    end
    if ~exist('FigPath','var') || isempty(FigPath)
        FigPath = pwd;
    end
    if ~exist('Name','var') || isempty(Name)
        Name = 'NetResults';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fig1 = figure;
    for node = 1:Net.NodeNum
        % 
        subplot(Net.NodeNum,2,(node-1)*2+1),plot(TS(node,1:Net.SF*2));
        set(gca,'xtick',0:Net.SF:Net.SF*2,'xticklabel',0:2)
        if node == Net.NodeNum, xlabel('time (s)'); end
        ylabel(['Node' num2str(node)],'fontweight','bold','fontsize',11);
        if node ==1, title('Temporal dynamic'); end

        subplot(Net.NodeNum,2,(node-1)*2+2);
        [Z,f] = pwelch(TS(node,:),Net.SF,[],[],Net.SF);
        plot(f,Z,'linewidth',2);
        if node==1, title('Power Spectrum Density');end
        if node == Net.NodeNum,  xlabel('Frequency(Hz)');end
        xlim([0 20]);
        set(gca,'yticklabel',[])
        %legend(['Signal #' num2str(Sig)])
    end
    if issave
        set(Fig1,'PaperPosition',[1 1 6 4.5])
        print(fullfile(FigPath,Name),'-r300','-dtiff')
    end
end

%%

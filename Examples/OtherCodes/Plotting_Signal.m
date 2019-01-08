clear; clc;
fs = 22;
%%
% SSVEP signal can be simulated using ModelSourceSignal with defined parameters, otherwise Roisignal function will generate a default two source SSVEP signal 
[outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[8 20],'signalHarmonic',{[2,0,0],[1,0,0]},'signalPhase',{[.1,0,0],[0,.3,0]});
FIG1 = figure;
[Z,f] = pwelch(mean(outSignal,2),SF,[],[],SF);
subplot(1,2,1),plot(mean(outSignal(1:100,:),2),'linewidth',2,'color','k');
xlabel('Time (S)','fontsize',fs);
set(gca,'fontsize',fs,'xtick',0:SF/2:SF,'xticklabel',0:.5:1,'yticklabel',[]);

subplot(1,2,2),plot(f,Z,'linewidth',2,'color','k');
xlim([0 30]);xlabel('Frequency (Hz)','fontsize',fs);
set(gca,'fontsize',fs,'yticklabel',[]);

set(FIG1,'paperposition',[1 1 8 2.8]);
print(fullfile('Figures','SignalModelSSEP.tif'),'-r300','-dtiff');
%%

load(fullfile('ExampleData','NetTimeSeries_3Nodes.mat'))
SF = 300;
FIG2 = figure;

[Z,f] = pwelch(TS_unconnect(3,1:2000),SF,[],[],SF);
subplot(1,2,1),plot(TS_unconnect(1,1:300),'linewidth',2,'color','k');
xlabel('Time (S)','fontsize',fs);
xlim([0 300])
set(gca,'fontsize',fs,'xtick',0:SF/2:SF,'xticklabel',0:.5:1,'yticklabel',[]);

subplot(1,2,2),plot(f,Z,'linewidth',2,'color','k');
xlim([0 30]);xlabel('Frequency (Hz)','fontsize',fs);
set(gca,'fontsize',fs,'yticklabel',[]);

set(FIG2,'paperposition',[1 1 8 2.8]);
print(fullfile('Figures','SignalModelARMA.tif'),'-r300','-dtiff');

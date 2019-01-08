clear; clc; close all;
AnatPath = 'C:\Users\Elhamkhanom\Documents\Codes\Git\mrC\Examples\ExampleData\anatomy';
%% plot two atlases
cmap1 = distinguishable_colors(25,[.5 .5 .5]);
CC1 = zeros(50,3);CC1(1:2:end,:)=cmap1;CC1(2:2:end,:)=cmap1;

cmap2 = distinguishable_colors(180,[.5 .5 .5]);
CC2 = zeros(360,3);CC2(1:2:end,:)=cmap2;CC2(2:2:end,:)=cmap2;
FIG = figure;
subplot(1,2,1),mrC.Simulate.VisualizeSourceRoi('nl-0048',AnatPath,'wang',[],'ventral','B',CC1,[]);
view (50,20)
subplot(1,2,2),mrC.Simulate.VisualizeSourceRoi('nl-0048',AnatPath,'glass',[],'ventral','B',CC2,[]);
view (50,20)
set(FIG,'PaperPosition',[1 1 8 3.5])
print('Figures/ATLASES.tif','-r300','-dtiff')

%% Plot V1, and signal modeling
FIG2 = figure;
CC1 = repmat([0 0 1],4,1);
S1 = subplot(1,2,1);
mrC.Simulate.VisualizeSourceRoi2('nl-0048',AnatPath,'wang',{'V1'},[],'B',CC1,'ventral');
view (50,20)
set(S1,'position',get(S1,'position')+[-0.1 -.1 .2 .15]);

load ExampleData/NetTimeSeries_3Nodes.mat;
subplot(2,2,2);
plot(TS_connect(1,1:300),'color','b','linewidth',2);
set(FIG2,'PaperPosition',[1 1 10 6]);
print('Figures/SeedSignal_nl0048.tif','-r300','-dtiff');


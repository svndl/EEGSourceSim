clear; clc;

mrCFolder = 'C:\Users\Elhamkhanom\Documents\Codes\Git\mrC';%'/Users/kohler/code/git/mrC';
addpath(genpath(mrCFolder));
FigPath = 'Figures'; % indicate where to save the resulting figures

%% Project of interest
DestPath = fullfile(mrCFolder,'Examples','ExampleData_Inverse');
AnatomyPath = fullfile(DestPath,'anatomy');
ProjectPath = fullfile(DestPath,'FwdProject');

%% modeling signal and noise and plot them
[outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[2 3],'signalHarmonic',{[2,0,1],[1,1,0]},'signalPhase',{[.1,0,.2],[0,.3,0]});

[RoiList,subIDs] = mrC.Simulate.GetRoiClass(ProjectPath,AnatomyPath);% 13 subjects with Wang atlab 
V1_RoiList = cellfun(@(x) {x.searchROIs('V1d','wang','L')},RoiList);
V3_RoiList = cellfun(@(x) {x.searchROIs('V3d','wang','L')},RoiList);
Seed_RoiList = cellfun(@(x,y,z) {mergROIs(x,y)},V1_RoiList,V3_RoiList);
%%
Noise.mu=2;
Noise.lambda = 1/length(outSignal)*2;

[EEGData1,EEGAxx1,SourceData,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',...
    AnatomyPath,'signalArray',outSignal,'signalsf',SF,'NoiseParams',Noise,'rois',...
    Seed_RoiList,'doSource' ,true,'subSelect','nl-0048','figFolder',FigPath,'VisualizeNoise',true);
    
%%
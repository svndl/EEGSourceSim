clear; clc;

SimFolder = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(SimFolder)));

if isempty(getpref('EEGSSim')) % if the Setup function hasn't been run yet
    EEGSourceSimSetUp();
elseif  isempty(getpref('EEGSSim','EEGSSimPath')) % if the EEGSourceSim hasn't been added to the path
    EEGSourceSimSetUp('EEGSourceSim',true);
end

%% Prepare the results folders
FigPath = fullfile(SimFolder,'Figures');
ResultPath = fullfile(SimFolder,'ResultData');
if ~exist(FigPath,'dir'),mkdir(FigPath);end
if ~exist(ResultPath,'dir'),mkdir(ResultPath);end

%% Prepare Project path and ROIs
if isempty(getpref('EEGSSim','AnatomyPath')) || isempty(getpref('EEGSSim','ProjectPath'))
    EEGSourceSimSetUp('EEGSourceSim',false,'Dataset',true);
end
AnatomyPath = getpref('EEGSSim','AnatomyPath');
ProjectPath = getpref('EEGSSim','ProjectPath');

%%
[RoiList,subIDs] = mrC.Simulate.GetRoiClass(ProjectPath,AnatomyPath);
Wangs = cellfun(@(x) {x.getAtlasROIs('wang')},RoiList);
Wangnums = cellfun(@(x) x.ROINum,Wangs)>0;


% define locations of sources
%--------------------------Cond1: V2d_R -----------------------------
Rois1 = cellfun(@(x) x.searchROIs('V2d','wang','R'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('LO1','wang','L'),RoiList,'UniformOutput',false);
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
do_new_data_generation = false;
% generate or read from disk
generated_date_filename = 'data_for_spatial_filter_test_2source_oneSubj.mat';

n_trials = 200;
Noise.lambda = 0 ; % noise only
[outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','ns',200,'signalFreq',[2 2],'harmonicAmps',{[2,0,1.5,0],[1,0, 1,0]},'harmonicPhases',{[0,0,0,0],[pi/2,0,pi/2,0]},'reliableAmps',[1,0],'nTrials',n_trials);

%%
for r = 1:10
    disp(['Run #' num2str(r)])
    T1 = clock;
    [EEGData_noise,EEGAxx_noise,EEGData_signal,EEGAxx_signal,~,masterList,subIDs,allSubjFwdMatrices,allSubjRois,Times] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,...
        'subSelect',subIDs,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',false,'cndNum',1,'nTrials',n_trials,'RedoMixingMatrices',true);
    T2 = clock;
    Times.overal = etime(T2,T1);
    Times.apha = sum(arrayfun(@(x) Times.noiset{x}.alpha,1:200));
    Times.pink = sum(arrayfun(@(x) Times.noiset{x}.pink,1:200));
    Times.white = sum(arrayfun(@(x) Times.noiset{x}.white,1:200));
    RunTimes2{r} = Times;
    clear Times;
end   

%%
for r = 1:10
    disp(['Run #' num2str(r)])
    T1 = clock;
    [EEGData_noise,EEGAxx_noise,EEGData_signal,EEGAxx_signal,~,masterList,subIDs,allSubjFwdMatrices,allSubjRois,Times] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,...
        'subSelect',subIDs,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',false,'cndNum',1,'nTrials',n_trials);%,'RedoMixingMatrices',true);
    T2 = clock;
    Times.overal = etime(T2,T1);
    Times.apha = sum(arrayfun(@(x) Times.noiset{x}.alpha,1:200));
    Times.pink = sum(arrayfun(@(x) Times.noiset{x}.pink,1:200));
    Times.white = sum(arrayfun(@(x) Times.noiset{x}.white,1:200));
    RunTimes{r} = Times;
    clear Times;
end    



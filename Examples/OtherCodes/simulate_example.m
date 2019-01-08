% This script generates an example simulated EEG with SSVEP signals
% for this script to run correctly you need three paths:
    % mrCpath: The latest mrC package
    % ProjectPath: Pointing to the mrC project folder 
    % AnatomyPath: pointing to the folder where anatomy data is 
                
% Elham Barzegaran 3/14/2018
% Latest Modification: 7/16/2018

%% Add latest mrC
clear;clc
mrCFolder = fileparts(fileparts(mfilename('fullpath')));%'/Users/kohler/code/git';
addpath(genpath(mrCFolder));

addpath(genpath('C:\Users\Elhamkhanom\Documents\Codes\Git\surfing'));% this tool can be found in github

%% mrC project and anatomy paths
% THIS PART OF SCRIPT IS FOR THE LAB COMPUTER, IF YOU ARE USING AN OUT OF STANFORD NETWORK
% COMPUTER, COMMENT THIS PART AND ONLY SET DestPath TO THE FOLDER YOU HAVE
% SAVED THE SAMPLE ANATOMY AND PROJECT FOLDER. YOU CAN DOWNLOAD THOSE FILES
% FROM MY DROPBOX: https://www.dropbox.com/sh/ypjrwjjw3003c2f/AABYO1JEjcdwkH3auBOon6UVa?dl=0
%
 DestPath = fullfile(mrCFolder,'Examples','ExampleData');
% 
% ProjectPath ='/Volumes/svndl/mrC_Projects/kohler/SYM_RT_LOCKED/SOURCE';
% % This is to make a portable copy of the project data (both anatomy and forward)
% [ProjectPath, AnatomyPath]  = mrC.Simulate.PrepareProjectSimulate(ProjectPath,DestPath ,'FwdFormat','mat');
% 
% ProjectPath = '/Volumes/svndl/mrC_Projects/kohler/SYM_16GR/SOURCE';
% % This is to make a portable copy of the project data (both anatomy and forward)
% mrC.Simulate.PrepareProjectSimulate(ProjectPath,DestPath ,'FwdFormat','mat');
% 
% 
% ProjectPath = '/Volumes/svndl/mrC_Projects/Att_disc_annulus/Source';
% % This is to make a portable copy of the project data (both anatomy and forward)
% mrC.Simulate.PrepareProjectSimulate(ProjectPath,DestPath ,'FwdFormat','mat');

%% Indicate project and anatomy folders and get ROIs for example subjects
% So far, 16 subject have Wang atlas ROIs in this dataset

AnatomyPath = fullfile(DestPath,'anatomy');
ProjectPath = fullfile(DestPath,'FwdProject');

% Pre-select ROIs
[RoiList,subIDs] = mrC.Simulate.GetRoiClass(ProjectPath,AnatomyPath);% 13 subjects with Wang atlab 
Wangs = cellfun(@(x) {x.getAtlasROIs('wang')},RoiList);
Wangnums = cellfun(@(x) x.ROINum,Wangs)>0;

%% SSVEP signal can be simulated using ModelSourceSignal with defined parameters, otherwise Roisignal function will generate a default two source SSVEP signal 
% a simple SSVEP signal...

[outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[2 3],'signalHarmonic',{[2,0,1],[1,1,0]},'signalPhase',{[.1,0,.2],[0,.3,0]});

%% simulating EEGs with different ROIs as different conditions
Noise.mu=2;
Noise.lambda = 1/length(outSignal);

%--------------------------Cond1: V2d_R, V3d_L-----------------------------
Rois1 = cellfun(@(x) x.searchROIs('V2d','wang','R'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('V3d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
[EEGData1,EEGAxx1,~,masterList1,subIDs1] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',1);
close all;

%--------------------------Cond2: V1d_L, V2d_L-----------------------------
Rois1 = cellfun(@(x) x.searchROIs('V1d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('V2d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
[EEGData2,EEGAxx2,~,masterList2,subIDs2] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',2);
close all;

%--------------------------Cond4: V1v_R, hV4_R-----------------------------
Rois1 = cellfun(@(x) x.searchROIs('V1v','wang','R'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('hV4','wang','R'),RoiList,'UniformOutput',false);% % wang ROI
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
[EEGData4,EEGAxx4,~,masterList4,subIDs4] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',4);
close all;

%--------------------------Cond3: V1d_L, V2d_L, LO1_L----------------------
[outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[2 3 5],'signalHarmonic',{[2,0,1],[1,1,0],[2, 1]},'signalPhase',{[.1,0,.2],[0,.3,0],[.5, .6]});
Rois1 = cellfun(@(x) x.searchROIs('V1d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('V2d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
Rois3 = cellfun(@(x) x.searchROIs('LO1','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
RoisI = cellfun(@(x,y) x.mergROIs(y),RoisI,Rois3,'UniformOutput',false);
[EEGData3,EEGAxx3,~,masterList3,subIDs3] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',3);
close all;

%% simulating EEGs with different ROIs as different conditions
Noise.mu=1;
[outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[2 3],'signalHarmonic',{[2,0,1],[1,1,0]},'signalPhase',{[.1,0,.2],[0,.3,0]});
Rois1 = cellfun(@(x) x.searchROIs('V2d','wang','R'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('V3d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);

%--------------------------Cond11: V2d_R, V3d_L,labda = 1/500--------------
Noise.lambda = 1/length(outSignal)*2;
[EEGData11,EEGAxx11,~,masterList11,subIDs11] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',11);
close all;
%--------------------------Cond12: V2d_R, V3d_L,labda = 1/750--------------
Noise.lambda = 1/length(outSignal)*(4/3);
[EEGData12,EEGAxx12,~,masterList12,subIDs12] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',12);
close all;
%--------------------------Cond13: V2d_R, V3d_L,labda = 1/1000--------------
Noise.lambda = 1/length(outSignal);
[EEGData13,EEGAxx13,~,masterList13,subIDs13] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',13);
close all;
%--------------------------Cond14: V2d_R, V3d_L,labda = 1/2000--------------
Noise.lambda = 1/length(outSignal)/2;
[EEGData14,EEGAxx14,~,masterList14,subIDs14] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',14);
close all;
%% visualize the results

Cond = 14; % condition to be visualized

fundfreq = [2 3 5];
eval(['EEGAxxP = EEGAxx' num2str(Cond) ';']);
eval(['masterListP = masterList' num2str(Cond) ';']);
sub1 = find(~cellfun(@isempty,EEGAxxP),1);
freq = 0:EEGAxxP{sub1}.dFHz:EEGAxxP{sub1}.dFHz*(EEGAxxP{sub1}.nFr-1); % frequncy labels, based on fft
for s = 1:length(EEGAxxP)
    if ~isempty(EEGAxxP{s})
        ASDEEG{s} = EEGAxxP{s}.Amp;% it is important which n is considered for fft
    end 
end

% Plot average over individuals
MASDEEG = mean(cat(4,ASDEEG{:}),4);
mrC.Simulate.PlotEEG(MASDEEG,freq,[],'average over all  ',masterListP,fundfreq(1:numel(masterListP)));








% This script modifies the script simulate_examples.m in order to test some
% new developments.
% for this script to run correctly you need three paths:
    % mrCpath: The latest mrC package
    % ProjectPath: Pointing to the mrC project folder 
    % AnatomyPath: pointing to the folder where anatomy data is 
                
% Elham Barzegaran 3/14/2018
% Modified by Sebastian Bosse: 8/2/2018

%% Add latest mrC
clear;clc
mrCFolder = fileparts(fileparts(mfilename('fullpath')));%'/Users/kohler/code/git';
%mrCFolder = '/Users/bosse/svn_dev/mrC_branched/mrC/' ;
addpath(genpath(mrCFolder));

%addpath(genpath('C:\Users\Elhamkhanom\Documents\Codes\Git\surfing'));% this tool can be found in github

%% mrC project and anatomy paths
% THIS PART OF SCRIPT IS FOR THE LAB COMPUTER, IF YOU ARE USING AN OUT OF STANFORD NETWORK
% COMPUTER, COMMENT THIS PART AND ONLY SET DestPath TO THE FOLDER YOU HAVE
% SAVED THE SAMPLE ANATOMY AND PROJECT FOLDER. YOU CAN DOWNLOAD THOSE FILES
% FROM MY DROPBOX: https://www.dropbox.com/sh/ypjrwjjw3003c2f/AABYO1JEjcdwkH3auBOon6UVa?dl=0
%
 %DestPath = fullfile(mrCFolder,'Examples','ExampleData');
 DestPath = '/export/data/eeg_simulation';
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
ProjectPath = fullfile(DestPath,'FwdProject2');

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
[EEGData1,EEGAxx1,~,masterList1,subIDs1] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',1,'nTrials',10);
close all;
%%
%--------------------------Cond2: V1d_L, V2d_L-----------------------------
Rois1 = cellfun(@(x) x.searchROIs('V1d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('V2d','wang','L'),RoiList,'UniformOutput',false);% % wang ROI
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
[EEGData2,EEGAxx2,~,masterList2,subIDs2] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',true,'cndNum',2,'nTrials',10);
close all;
%%


%% 
do_whitening = true;
freq = EEGAxx1{1}.dFHz*[0:(EEGAxx1{1}.nFr-1)];

sigFreqs = {[2,6],[3,6]} ;
fundFreqs = [2,3] ;
masterListP = [masterList1, masterList2] ;
%%
MASDEEG_raw_comp = mean(EEGAxx1{1}.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_raw_comp,freq,[],'raw  condition 1',masterList1,fundFreqs);
%%
MASDEEG_raw_comp = mean(EEGAxx2{1}.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_raw_comp,freq,[],'raw  condition 2',masterList2,fundFreqs);
%% PCA on full freq range
[PCAAxx,W,A,D] = mrC.SpatialFilters.PCA(EEGAxx1{1},'freq_range',freq);
MASDEEG_pca_comp = mean(PCAAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_pca_comp,freq,[],'PCA on full freq. range ',masterList1,fundFreqs,[],[],[],[],[],A,W);

%% PCA on 1st stim freq
[PCAAxx,W,A,D] = mrC.SpatialFilters.PCA(EEGAxx1{1},'freq_range',fundFreqs(1));
MASDEEG_pca_comp = mean(PCAAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_pca_comp,freq,[], sprintf('PCA on %0.1f Hz', fundFreqs(1)), masterList1,fundFreqs,[],[],[],[],[],A,W);

%% PCA on 2st stim freq
[PCAAxx,W,A,D] = mrC.SpatialFilters.PCA(EEGAxx1{1},'freq_range',fundFreqs(2));
MASDEEG_pca_comp = mean(PCAAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_pca_comp,freq,[], sprintf('PCA on %0.1f Hz', fundFreqs(2)), masterList1,fundFreqs,[],[],[],[],[],A,W);

%% SSD on 1st stim freq
[SSDAxx,W,A,D] = mrC.SpatialFilters.SSD(EEGAxx1{1},fundFreqs(1),'do_whitening',true);
MASDEEG_ssd_comp = mean(SSDAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_ssd_comp,freq,[],sprintf('SSD on %0.1f Hz', fundFreqs(1)),masterList1,fundFreqs,[],[],[],[],[],A,W);

%% SSD on 2nd stim freq
[SSDAxx,W,A,D] = mrC.SpatialFilters.SSD(EEGAxx1{1},fundFreqs(2),'do_whitening',true);
MASDEEG_ssd_comp = mean(SSDAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_ssd_comp,freq,[],sprintf('SSD on %0.1f Hz', fundFreqs(2)),masterList1,fundFreqs,[],[],[],[],[],A,W);

%% CSP
[CSPAxx,W,A,D] = mrC.SpatialFilters.CSP({EEGAxx1{1},EEGAxx2{1}},'freq_range',[2,3],'do_whitening',do_whitening);
class_vis = 1;
MASDEEG_csp_comp = mean(CSPAxx{class_vis}.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_csp_comp,freq,[],'CSP on fundFreqs ',masterList1,fundFreqs,[],[],[],[],[],A,W);

%% CSP
[CSPAxx,W,A,D] = mrC.SpatialFilters.CSP({EEGAxx1{1},EEGAxx2{1}},'do_whitening',do_whitening);
class_vis = 1;
MASDEEG_csp_comp = mean(CSPAxx{class_vis}.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_csp_comp,freq,[],'CSP on all freqs ',masterList1,fundFreqs,[],[],[],[],[],A,W);

%% RCA
[RCAAxx,W,A,D] = mrC.SpatialFilters.RCA(EEGAxx1{1},'freq_range',[2,3],'do_whitening',do_whitening);
MASDEEG_rca_comp = mean(RCAAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_rca_comp,freq,[],'RCA on fundFreqs ',masterList1,fundFreqs,[],[],[],[],[],A,W);

%% RCA
[RCAAxx,W,A,D] = mrC.SpatialFilters.RCA(EEGAxx1{1},'freq_range',[2,3],'do_whitening',do_whitening, 'model_type','cartesian');
MASDEEG_rca_comp = mean(RCAAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_rca_comp,freq,[],'RCA on fundFreqs, indep. cplx comps. ',masterList1,fundFreqs,[],[],[],[],[],A,W);

%% RCA
[RCAAxx,W,A,D] = mrC.SpatialFilters.RCA(EEGAxx1{1},'do_whitening',do_whitening);
MASDEEG_rca_comp = mean(RCAAxx.Amp,3);
mrC.Simulate.PlotEEG(MASDEEG_rca_comp,freq,[],'RCA on all freqs ',masterList1,fundFreqs,[],[],[],[],[],A,W);

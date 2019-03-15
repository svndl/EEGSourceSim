% This Script generates simulation for spatial filters
% Author: Sebastian Bosse,
% Modfied: Elham Barzegaran, 1/2019


%% Add latest mrC
clear;clc
SimFolder = fileparts(pwd);
addpath(genpath(SimFolder));

%% To be modified later
if false % SBs setup
    addpath('../External/tools/BrewerMap/')
    %
    DataPath = '/export/data/';
    DestPath = fullfile(DataPath,'eeg_simulation');
    AnatomyPath = fullfile(DestPath,'anatomy');
    ProjectPath = fullfile(DestPath,'FwdProject');
else
    DestPath = fullfile(SimFolder,'Examples','ExampleData_Inverse');
    AnatomyPath = fullfile(DestPath,'anatomy');
    ProjectPath = fullfile(DestPath,'FwdProject');
end

%% Prepare the results folders
FigPath = fullfile(SimFolder,'Examples','Figures');
ResultPath = fullfile(SimFolder,'Examples','ResultData');
if ~exist(fullfile(FigPath),'dir'),mkdir(FigPath);end
if ~exist(fullfile(ResultPath),'dir'),mkdir(ResultPath);end

%%
% Pre-select ROIs
[RoiList,subIDs] = mrC.Simulate.GetRoiClass(ProjectPath,AnatomyPath);% 13 subjects with Wang atlab 
Wangs = cellfun(@(x) {x.getAtlasROIs('wang')},RoiList);
Wangnums = cellfun(@(x) x.ROINum,Wangs)>0;

% define noise properties
Noise.mu.pink=2;
Noise.mu.alpha=2;
Noise.mu.sensor=2;
%Noise.distanceType = 'geodesic';

% define locations of sources
%--------------------------Cond1: V2d_R -----------------------------
Rois1 = cellfun(@(x) x.searchROIs('V2d','wang','R'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('LO1','wang','L'),RoiList,'UniformOutput',false);
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
do_new_data_generation = false;
% generate or read from disk
generated_date_filename = 'data_for_spatial_filter_test2_2source_allSubj.mat';
%generated_date_filename = 'data_for_spatial_filter_test_2source_all_subjects.mat';
%if ~exist('data_for_spatial_filter_test_2source.mat','file') || do_new_data_generation
if ~exist(generated_date_filename,'file') || do_new_data_generation
    n_trials = 20;
    Noise.lambda = 0 ; % noise only
    [outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','ns',200,'signalFreq',[2 2],'harmonicAmps',{[2,0,1.5,0],[1,0, 1,0]},'harmonicPhases',{[0,0,0,0],[pi/2,0,pi/2,0]},'reliableAmps',[1,0],'nTrials',n_trials,'reliableAmp',[1 0]);
    [EEGData_noise,EEGAxx_noise,EEGData_signal,EEGAxx_signal,~,masterList,subIDs,allSubjFwdMatrices,allSubjRois] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',false,'cndNum',1,'nTrials',n_trials);%,'RedoMixingMatrices',true);
    save(fullfile(ResultPath,generated_date_filename),'-v7.3');
else
    load(fullfile(ResultPath,generated_date_filename))
end
%%
% clean empty EEG data
idxs_empty = cellfun(@isempty,EEGAxx_noise) ;
EEGData_noise(idxs_empty) = [] ;
EEGData_signal(idxs_empty) = [] ;
EEGAxx_noise(idxs_empty) = [] ;
EEGAxx_signal(idxs_empty) = [] ;
subIDs(idxs_empty) = [] ;

%%
% mix signal and nose according to SNR and  convert to Axx
opt.signalFF=FundFreq(1) ;
opt.signalsf=SF ;
opt.cndNum = 1;

% note that this SNR is defined over the full spectrum, while the signal is
% narrowbanded
% SNR parameters
F1 = EEGAxx_signal{1,1}.i1F1;

Db_list = [-20:5:10] ; 
Lambda_list = 10.^(Db_list/10) ;

% spatial filter test parameters
fund_freq_idx = 1 ;            
nTrials = 20;
nDraws = 40 ;
n_comps = 3 ;
thisFundFreq = FundFreq(fund_freq_idx) ;

subs = num2cell(1:min(sum(not(idxs_empty)),10)) ; %%%%%% SUBJECTS TO SELECT
subNames = cellfun(@num2str,subs,'uni',false);

%EEGData_noise = cellfun(@(x) x(:,:,1:200),EEGData_noise,'uni',false); % reduce data size
%EEGAxx_noise = cellfun(@(x) x.SelectTrials(1:200),EEGAxx_noise,'uni',false);

rois = allSubjRois(not(idxs_empty)) ;
rois = rois(cell2mat(subs)) ;
fwdMatrix = allSubjFwdMatrices(not(idxs_empty)) ;
fwdMatrix = fwdMatrix(cell2mat(subs)) ;
Source_pattern = zeros(size(fwdMatrix{1},1),rois{1}.ROINum,numel(subs)) ;
for sub = 1:numel(subs)
    for roi_idx = 1:rois{1}.ROINum
        % assuming uniform activation
        Source_pattern(:,roi_idx,sub) = sum(fwdMatrix{subs{sub}}(:,rois{subs{sub}}.ROIList(roi_idx).meshIndices ),2)  ;
    end
end

% narrow-band normalization
% harmonics considered for normalization
power_norm_harmonics = [1,2,3,4] ;

f = [0:size(EEGData_signal{1},1)-1] *SF/size(EEGData_signal{1},1) ;

[~,power_norm_signal_freq_idxs]=intersect(f,power_norm_harmonics*thisFundFreq) ;
[~,power_norm_noise_freq_idxs]=intersect(f,power_norm_harmonics*thisFundFreq) ;

for subj_idx = 1:length(subIDs)
    %
    spec_noise = fft(EEGData_noise{subj_idx},[],1);
    spec_signal = fft(EEGData_signal{subj_idx},[],1);
    power_noise = mean(mean(abs(spec_noise(power_norm_noise_freq_idxs,:,:)).^2)) ; % mean noise power per trial
    power_signal= mean(mean(abs(spec_signal(power_norm_signal_freq_idxs,:,:)).^2)) ; % mean noise power per trial
    EEGData_noise{subj_idx} = EEGData_noise{subj_idx}./sqrt(power_noise);
    EEGData_signal{subj_idx} = EEGData_signal{subj_idx}./sqrt(power_signal);
end

% free some memory
clear spec_noise spec_signal fwdMatrix RoiList allSubjFwdMatrices allSubjRois Rois1 Rois2 Wangs Wangnums

if size(EEGData_signal{subj_idx},3)>= nDraws*nTrials % enough trials for disjoint draws
    trial_idxs = 1:nDraws*nTrials ;
else
    trial_idxs = randi( size(EEGData_signal{subj_idx},3), 1,nDraws*nTrials) ; % not enough trials
end
trial_idxs_per_draw = reshape(trial_idxs,nDraws,[]) ;

decomp_methods = {'pca','ssd','csp','rca'} ;
for decomp_method_idx = 1:length(decomp_methods)
    this_decomp_method = decomp_methods{decomp_method_idx};
    for s = 1:numel(subs)
        ang_errs.(this_decomp_method){s} = zeros(2,128,numel(Lambda_list),nDraws) ;
        residuals.(this_decomp_method){s} = zeros(2,128,numel(Lambda_list),nDraws) ;
    end
end
    
ref_signals = squeeze(mean(outSignal(1:100,:,:),2)) ;
%ref_signals = squeeze(outSignal(1:100,:,:));

for nLambda_idx = 1:numel(Lambda_list)
    lambda = Lambda_list(nLambda_idx);
    disp(['Generating EEG by adding signal and noise: SNR = ' num2str(lambda)]);
    for subj_idx = 1:length(subIDs)
        EEGData{subj_idx} = sqrt(lambda/(1+lambda))*EEGData_signal{subj_idx} + sqrt(1/(1+lambda)) * EEGData_noise{subj_idx} ;
        EEGAxx{subj_idx} = mrC.Simulate.CreateAxx(EEGData{subj_idx},opt) ;
    end

    % test spatial filters  
    for s = 1:numel(subs)
        display(['Calculating spatial filters for subject:' subNames{s}]);
        % TODO: what is this?
        subj_idx = subs{s};
        
        source_pattern = Source_pattern(:,:,subj_idx );
        considered_harms = power_norm_harmonics ;

        for draw_idx = 1:nDraws
            trial_idxs = trial_idxs_per_draw(draw_idx,:);

            thisAxx = cellfun(@(x) x.SelectTrials(trial_idxs),EEGAxx(subj_idx),'uni',false);
            T = thisAxx{1};
            for sub = 2:numel(thisAxx)
                T = T.MergeTrials(thisAxx{sub});
            end
            thisAxx = T;
            thisTempMean = mean(mean(thisAxx.Wave,3),1) ;

            % make sure the no-stimulation condition does not see the same
            % noise component
            noise_trial_idxs = trial_idxs_per_draw(mod(draw_idx,nDraws)+1,:);
            
            thisNoiseAxx = cellfun(@(x) x.SelectTrials(trial_idxs),EEGAxx_noise(subj_idx),'uni',false);
            T = thisNoiseAxx{1};
            for sub = 2:numel(thisNoiseAxx)
                T = T.MergeTrials(thisNoiseAxx{sub});
            end
            thisNoiseAxx = T;

            for decomp_method_idx = 1:length(decomp_methods)
                this_decomp_method = decomp_methods{decomp_method_idx};

                if strcmpi(this_decomp_method,'pca')
                    [thisDecompAxx,thisW,thisA,thisD] = mrC.SpatialFilters.PCA(thisAxx,'freq_range',thisFundFreq*considered_harms);
                elseif strcmpi(this_decomp_method,'pca_cart')
                    [thisDecompAxx,thisW,thisA,thisD] = mrC.SpatialFilters.PCA(thisAxx,'freq_range',thisFundFreq*considered_harms,'model_type','cartesian');
                elseif strcmpi(this_decomp_method,'fullfreqPca')
                    [thisDecompAxx,thisW,thisA,thisD] = mrC.SpatialFilters.PCA(thisAxx);
                elseif strcmpi(this_decomp_method,'fullfreqPca')
                    [thisDecompAxx,thisW,thisA,thisD] = mrC.SpatialFilters.PCA(thisAxx,'freq_range',[1:50]);
                elseif strcmpi(this_decomp_method,'tpca')
                    [thisDecompAxx,thisW,thisA,thisD] = mrC.SpatialFilters.tPCA(thisAxx);
                elseif strcmpi(this_decomp_method,'ssd')
                    [thisDecompAxx,thisW,thisA,thisD]= mrC.SpatialFilters.SSD(thisAxx,thisFundFreq*considered_harms,'do_whitening',true);
                elseif strcmpi(this_decomp_method,'rca')
                    [thisDecompAxx,thisW,thisA,thisD] = mrC.SpatialFilters.RCA(thisAxx,'freq_range',thisFundFreq*considered_harms,'do_whitening',true);

                elseif strcmpi(this_decomp_method,'csp')
                    [theseDecompAxxs,thisW,thisA,thisD] = mrC.SpatialFilters.CSP({thisAxx,thisNoiseAxx},'freq_range',thisFundFreq*considered_harms,'do_whitening',true);
                    thisDecompAxx=theseDecompAxxs{1};
                end
                
                % reduce size of decomposition result to save memory,
                % otherwise matlab seems to crash
                
                thisW = thisW(:,1:5) ;
                thisA = thisA(:,1:5) ;
                thisD = thisD(1:5) ;
                thisDecompAxx.Amp = thisDecompAxx.Cos(:,1:5,:) ;
                thisDecompAxx.Cos = thisDecompAxx.Cos(:,1:5,:) ;
                thisDecompAxx.Sin = thisDecompAxx.Sin(:,1:5,:) ;
                thisDecompAxx.Wave = thisDecompAxx.Wave(:,1:5,:) ;
                
                
                % calcualte angular error and turn activations, filters and
                % signals if approriate

                these_ang_errs = 180/pi* acos((source_pattern'*thisA)./sqrt(repmat(sum(source_pattern.^2)',[1 size(thisA,2)]).*repmat(sum(thisA.^2),[size(source_pattern,2) 1]))) ;
                for cand_idx =  find(these_ang_errs>90)' % candidates for turning topopgraphies
                    [source_idx,comp_idx] =ind2sub(size(these_ang_errs), cand_idx) ;
                    if (180-these_ang_errs(source_idx,comp_idx)) <= min([these_ang_errs(:,comp_idx);180-these_ang_errs(:,comp_idx)])
 % turn activation, filters and signal if turned angular error is the smallest
                        these_ang_errs(source_idx,comp_idx) = 180-these_ang_errs(source_idx,comp_idx) ;
                        thisW(:,comp_idx) = -1*thisW(:,comp_idx) ;
                        thisA(:,comp_idx) = -1*thisA(:,comp_idx) ;
                        thisDecompAxx.Cos(:,comp_idx,:) = -1*thisDecompAxx.Cos(:,comp_idx,:) ;
                        thisDecompAxx.Sin(:,comp_idx,:) = -1*thisDecompAxx.Sin(:,comp_idx,:) ;
                        thisDecompAxx.Wave(:,comp_idx,:) = -1*thisDecompAxx.Wave(:,comp_idx,:) ;
                    end
                end
                
                Axx_compspace.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisDecompAxx ;
                W.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisW ;
                A.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisA ;
                D.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisD ;
                
                
                ang_errs.(this_decomp_method){s}(:,1:size(these_ang_errs,2),nLambda_idx,draw_idx) = these_ang_errs ;
                

                %calculate snrs assuming ssveps, mean over all trials
                freqs = [0:thisDecompAxx.nFr]*thisDecompAxx.dFHz;
                signal_freq_idxs = find(ismember(freqs,thisFundFreq*considered_harms));
                noise_freq_idxs = reshape([signal_freq_idxs-1;signal_freq_idxs+1],1,[]) ;

                snrs.(this_decomp_method){s}(1:size(thisA,2),nLambda_idx,draw_idx)=mean(mean(thisDecompAxx.Amp(signal_freq_idxs,:,:).^2)./mean(thisDecompAxx.Amp(noise_freq_idxs,:,:).^2),3);
                  % residual (mse over sampmles averaged and trials)
                % calculate residuals as mse over samples and trials
                est_signal = thisDecompAxx.Wave ;               

                corr_per_trial = zeros(2,thisDecompAxx.nCh,thisDecompAxx.nTrl) ;
                for trial_idx = 1:thisDecompAxx.nTrl
                    corr_per_trial(:,:,trial_idx) = corr(ref_signals, est_signal(:,:,trial_idx)) ;
                end
                residuals.(this_decomp_method){s}(:,1:thisDecompAxx.nCh,nLambda_idx,draw_idx) = mean(corr_per_trial.^2,3); % average residual over trials
                
            end
            snrs_orig{s}(:,nLambda_idx,draw_idx)=mean(mean(thisAxx.Amp(signal_freq_idxs,:,:).^2)./mean(thisAxx.Amp(noise_freq_idxs,:,:).^2),3);

        end
    end

    EEGData = {};
    EEGAxx = {} ;
end
%%
do_save_plots = true;

%% scalp plots
% settings for data
Subject_idx = 1;
source_pattern = Source_pattern(:,:,Subject_idx);
n_comps = 2;
decomp_methods = {'pca','rca','ssd'};

% settings for plot
fig_width = 8.5 ;
top_margin = 0.5 ;
left_margin = 0 ;
text_height = 0.25 ;
text_width = 1 ;
hor_comp_dist = 0.2 ;

sbpl_width = (fig_width-text_width-hor_comp_dist*length(decomp_methods))/(1+n_comps*length(decomp_methods));
sbpl_height = sbpl_width;
fig_height = top_margin+2*text_height + length(Lambda_list)* sbpl_height ;

fig_scalp_plots = figure();
set(fig_scalp_plots,'Units','Inches','position',[10, 10, fig_width, fig_height],'color','w');

for source_idx = 1:2
    x = left_margin ;
    y = fig_height - (top_margin+text_height + 0.5*length(Lambda_list)*sbpl_height + (source_idx-1)* (sbpl_height+text_height) ) ;
    ax = axes('parent',fig_scalp_plots,    'Units','Inches','Position',[x,y,sbpl_width,sbpl_height]) ;

    Topo = source_pattern(:,source_idx);
    %if abs(min(Topo))>max(Topo), Topo = -1*Topo;end
    
    this_title = sprintf('%s %i','Source ',source_idx) ;
    mrC.Simulate.PlotScalp(Topo,this_title);
end

for this_decomp_method_idx = 1:length(decomp_methods)
    this_decomp_method = decomp_methods{this_decomp_method_idx};
    for this_comp_idx = 1:2
        x = left_margin+hor_comp_dist+   this_comp_idx* sbpl_width +(this_decomp_method_idx-1)*(2*sbpl_width+hor_comp_dist) ;

        for nLambda_idx = 1:length(Lambda_list)
            y = fig_height - (top_margin+text_height + nLambda_idx*sbpl_height) ;
            ax = axes('parent',fig_scalp_plots,    'Units','Inches','Position',[x,y,sbpl_width,sbpl_height]) ;
            Topo = A.(this_decomp_method){Subject_idx}{nLambda_idx}{1}(:,this_comp_idx);
            if nLambda_idx == 1
                this_title = sprintf('%s %i','Comp ',this_comp_idx) ;
                mrC.Simulate.PlotScalp(Topo,this_title);
                if this_comp_idx == 1
                    ax = axes('parent',fig_scalp_plots,    'Units','Inches','Position',[x+sbpl_width,y+sbpl_height+text_height,0,0]) ;
                    text(0,0,upper(this_decomp_method),'Units','Inches', 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12,'FontWeight', 'bold');axis off
                end
            else
                mrC.Simulate.PlotScalp(Topo);
            end
            if (this_decomp_method_idx==length(decomp_methods)) && (this_comp_idx==2)
                ax = axes('parent',fig_scalp_plots,    'Units','Inches','Position',[x+sbpl_width+hor_comp_dist,y+0.5*sbpl_height,0,0]) ;
                this_text = sprintf('%i dB',Db_list(nLambda_idx)) ;
                text(0,0,this_text,'Units','Inches', 'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12,'FontWeight', 'bold');axis off
            end
            
        end
    end
end

if do_save_plots
    export_fig(fig_scalp_plots,fullfile(FigPath,'SpatialFilters_topographies'),'-pdf','-nocrop');
    close(fig_scalp_plots);
end

%
%% second figure
%  plots angular error of topographies

FS = 12;
fig_width = 8.5 ;
fig_height = fig_width*5/8/2 ;
fig_err_all = figure ;


colors = brewermap(4,'Set2') ;
markers = {'-o','-o'};
legend_entries = upper(decomp_methods);
for comp_idx = 1:2 
    subplot(3,2,comp_idx)
    for decomp_method_idx=1:length(decomp_methods)
        this_decomp_method = decomp_methods{decomp_method_idx};
        combined_err_angles.(this_decomp_method) = cat(5, ang_errs.(this_decomp_method){:});
        this_avg_err_angles = squeeze(mean(mean(combined_err_angles.(this_decomp_method)(comp_idx,comp_idx,:,:,:),4),5)) ;
        plot(Db_list,this_avg_err_angles ,markers{comp_idx},'LineWidth',2,'MarkerSize',10,'color',colors(decomp_method_idx,:));
        hold on
    end
    title(sprintf('Comp %i ( Relative to source %i )', comp_idx, comp_idx),'FontSize',FS)
    if comp_idx==1
        ylabel('Angular Error [Degrees]','FontSize',FS);
    else
        set(gca,'position', get(gca,'position')+[-.03 0 0 0]);
    end
    %xlabel('SNR [dB]','fontsize',FS);
    set(gca,'FontSize',FS,'xtick',-20:5:10)
end
[~, hobj, ~, ~] = legend(legend_entries,'Location','southwest');


for comp_idx = 1:2
    subplot(3,2,comp_idx+2)
    for decomp_method_idx=1:length(decomp_methods)
        this_decomp_method = decomp_methods{decomp_method_idx};
        combined_residuals.(this_decomp_method) = cat(5, residuals.(this_decomp_method){:}) ;
        this_avg_residuals = squeeze(mean(mean(combined_residuals.(this_decomp_method)(comp_idx,comp_idx,:,:,:),4),5)) ;
        plot(Db_list,this_avg_residuals ,markers{comp_idx},'LineWidth',2,'MarkerSize',10,'color',colors(decomp_method_idx,:));
        hold on
    end
    title(sprintf('Comp %i ( Relative to source %i )', comp_idx, comp_idx),'FontSize',FS)
    ylim([0,1])
    if comp_idx == 1
        ylabel('Normalized residual [1-R^2]','FontSize',FS);
    else
        set(gca,'position', get(gca,'position')+[-.03 0 0 0]);
    end
    %xlabel('SNR [dB]','fontsize',FS);
    set(gca,'FontSize',FS,'xtick',-20:5:10)
end

for comp_idx = 1:2
    subplot(3,2,comp_idx+4)
    for decomp_method_idx=1:length(decomp_methods)
        this_decomp_method = decomp_methods{decomp_method_idx};
        combined_snrs.(this_decomp_method)  =cat(4,snrs.(this_decomp_method){:}) ;
        %this_avg_snrs = mean(mean(10*log10(combined_snrs.(this_decomp_method)),3),4);
        this_avg_snrs = squeeze(mean(mean(10*log10(combined_snrs.(this_decomp_method)(comp_idx,:,:,:)),3),4));
        plot(Db_list,this_avg_snrs ,markers{comp_idx},'LineWidth',2,'MarkerSize',10,'color',colors(decomp_method_idx,:));
        hold on
    end
    title(sprintf('Comp %i', comp_idx),'FontSize',FS)
    xlabel('SNR [dB]','fontsize',FS);
    if comp_idx == 1
        ylabel('Reconstruction SNR [dB]','FontSize',FS); 
    else
        set(gca,'position', get(gca,'position')+[-.03 0 0 0]);
    end
    set(gca,'FontSize',FS,'xtick',-20:5:10)
    ylim([-5 35])
end


set(fig_err_all,'Units','Inches','position',[10, 10, 9, 9],'color','w');
set(fig_err_all,'paperposition',[10, 10, 9, 9]);

if do_save_plots
    print(fullfile(FigPath,'Spatial_Filters_ErrorPlots'),'-r300','-dtiff');
    export_fig(fig_err_all,fullfile(FigPath,'Spatial_Filters_ErrorPlots'),'-pdf');
    close(fig_err_all);
end














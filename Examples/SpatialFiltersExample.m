% This Script generates simulation for spatial filters
% Author: Sebastian Bosse,
%Modfied: Elham Barzegaran, 1/2019


%% Add latest mrC
clear;clc
SimFolder = fileparts(pwd);
addpath(genpath(SimFolder));

%% To be modified later
if false % SBs setup
    addpath('../tools/BrewerMap/')
    %
    DataPath = '/export/data/';
    DestPath = fullfile(DataPath,'eeg_simulation');
    AnatomyPath = fullfile(DestPath,'anatomy');
    ProjectPath = fullfile(DestPath,'FwdProject2');
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

% define locations of sources
%--------------------------Cond1: V2d_R -----------------------------
Rois1 = cellfun(@(x) x.searchROIs('V2d','wang','R'),RoiList,'UniformOutput',false);% % wang ROI
Rois2 = cellfun(@(x) x.searchROIs('LO1','wang','L'),RoiList,'UniformOutput',false);
RoisI = cellfun(@(x,y) x.mergROIs(y),Rois1,Rois2,'UniformOutput',false);
do_new_data_generation = false;
% generate or read from disk
if ~exist('data_for_spatial_filter_test_2source.mat','file') || do_new_data_generation
    n_trials = 5;
    Noise.lambda = 0 ; % noise only
    [outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[2 2],'signalHarmonic',{[2,0,1.5,0],[1.5,0, 2,0]},'signalPhase',{[.1,0,.1,0],[pi/2+.1,0,pi/2+.1,0]});
    [EEGData_noise,EEGAxx_noise,EEGData_signal,EEGAxx_signal,~,masterList,subIDs,allSubjFwdMatrices,allSubjRois] = mrC.Simulate.SimulateProject(ProjectPath,'anatomyPath',AnatomyPath,'signalArray',outSignal,'signalFF',FundFreq,'signalsf',SF,'NoiseParams',Noise,'rois',RoisI,'Save',false,'cndNum',1,'nTrials',n_trials);
    save(fullfile(ResultPath,'data_for_spatial_filter_test2_2source.mat'));
    save(fullfile(ResultPath,'data_for_spatial_filter_test_2source.mat'),'EEGAxx_noise','EEGData_noise','-v7.3');
else
    load(fullfile(ResultPath,'data_for_spatial_filter_test2_2source.mat'))
    load(fullfile(ResultPath,'data_for_spatial_filter_test_2source.mat'))
end

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
numTrials_list = 20;%2.^[1,2,3,4,5,6,7];     
nDraws = 20 ;
n_comps = 3 ;
thisFundFreq = FundFreq(fund_freq_idx) ;

subs = num2cell(1:min(length(subIDs),10)) ; %%%%%% SUBJECTS TO SELECT
subNames = cellfun(@num2str,subs(1:min(length(subIDs),10)),'uni',false);

EEGData_noise = cellfun(@(x) x(:,:,1:200),EEGData_noise,'uni',false); % reduce data size
EEGAxx_noise = cellfun(@(x) x.SelectTrials(1:200),EEGAxx_noise,'uni',false);

rois = allSubjRois(cell2mat(subs)) ;
fwdMatrix = allSubjFwdMatrices(cell2mat(subs)) ;
Source_pattern = zeros(size(fwdMatrix{1},1),rois{1}.ROINum,numel(subs)) ;
for sub = 1:numel(subs)
    for roi_idx = 1:rois{1}.ROINum
        % assuming uniform activation
        Source_pattern(:,roi_idx,sub) = sum(fwdMatrix{subs{sub}}(:,rois{subs{sub}}.ROIList(roi_idx).meshIndices ),2)  ;
    end
end

% narrow-band normalization
% harmonics considered for normalization
snr_harmonics = [1,2,3,4] ;

f = [0:length(EEGData_signal{1})-1] *SF/length(EEGData_signal{1}) ;
[~,snr_freq_idxs]=intersect(f,snr_harmonics*thisFundFreq) ;

for subj_idx = 1:length(subIDs)
    %
    spec_noise = fft(EEGData_noise{subj_idx},[],1);
    spec_signal = fft(EEGData_signal{subj_idx},[],1);
    power_noise = mean(mean(abs(spec_noise(snr_freq_idxs,:,:)).^2)) ; % mean noise power per trial
    power_signal= mean(mean(abs(spec_signal(snr_freq_idxs,:,:)).^2)) ; 
    EEGData_noise{subj_idx} = EEGData_noise{subj_idx}./sqrt(power_noise);
    EEGData_signal{subj_idx} = EEGData_signal{subj_idx}./sqrt(power_signal);
end


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
        subj_idx = subs{s};
        
        source_pattern = Source_pattern(:,:,subj_idx );

        decomp_methods = {'pca','ssd','csp','rca'} ;
        considered_harms = snr_harmonics ;

        for nTrial_idx = 1:length(numTrials_list)
            nUsedTrials = numTrials_list(nTrial_idx);
            for draw_idx = 1:nDraws
                fprintf('nUsedTrials = %i, draw_idx = %i \n',nUsedTrials, draw_idx) ;
                random_numbers = randperm(EEGAxx_noise{1}.nTrl) ;
                % todo: improve indexing to avoid identical trials in different
                % draws. requires simulation of enough trials
                % for starters: just take random trials
                trial_idxs = random_numbers(1:nUsedTrials);

                thisAxx = cellfun(@(x) x.SelectTrials(trial_idxs),EEGAxx(subj_idx),'uni',false);
                T = thisAxx{1};
                for sub = 2:numel(thisAxx)
                    T = T.MergeTrials(thisAxx{sub});
                end
                thisAxx = T;
                thisTempMean = mean(mean(thisAxx.Wave,3),1) ;

                % make sure the no-stimulation condition does not see the same
                % noise component
                noise_trial_idxs = random_numbers(nUsedTrials+1:2*nUsedTrials);
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
                    for i = 1:size(thisA,2)
                        if source_pattern(:,1)'*thisA(:,i)<0
                            thisA(:,i) = thisA(:,i)*-1 ;
                        end
                    end
                    Axx_compspace.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisDecompAxx ;
                    W.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisW ;
                    A.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisA ;
                    D.(this_decomp_method){s}{nLambda_idx}{draw_idx} = thisD ;

                    % metrics for first 2 components
                    % calculate error angles
                    freqs = [0:thisDecompAxx.nFr]*thisDecompAxx.dFHz;
                    signal_freq_idxs = find(ismember(freqs,thisFundFreq*considered_harms));
                    noise_freq_idxs = reshape([signal_freq_idxs-1;signal_freq_idxs+1],1,[]) ;

                    %err_angles.(this_decomp_method)(comp_idx,nTrial_idx,draw_idx) = 180/pi* acos(abs(source_pattern(:,1)'*thisA(:,comp_idx))/sqrt(sum(source_pattern(:,1).^2)*sum(thisA(:,comp_idx).^2))) ;
                    ERAngle = 180/pi* acos(abs(source_pattern'*thisA)./sqrt(repmat(sum(source_pattern.^2)',[1 size(thisA,2)]).*repmat(sum(thisA.^2),[size(source_pattern,2) 1]))) ;
                    ERAngle = ERAngle(1:2,1:2);
                    if ERAngle (1,1)> ERAngle (1,2)% if the order of components is flipped
                            ERAngle = ERAngle(1:2,2:-1:1);
                    end
                    err_angles.(this_decomp_method){s}(:,1:2,nLambda_idx,draw_idx) = ERAngle;
                    err_angles.(this_decomp_method){s}(:,:,nLambda_idx,draw_idx)=...
                            real(err_angles.(this_decomp_method){s}(:,:,nLambda_idx,draw_idx));


                    %calculate snrs assuming ssveps, mean over all trials
                    snrs.(this_decomp_method){s}(1:size(thisA,2),nLambda_idx,draw_idx)=mean(mean(thisDecompAxx.Amp(signal_freq_idxs,:,:).^2)./mean(thisDecompAxx.Amp(noise_freq_idxs,:,:).^2),3);
                      % residual (mse over sampmles averaged and trials)
                    % calculate residuals as mse over samples and trials
                    est_signal = thisDecompAxx.Wave ;               
                    ref_signal = outSignal(1:100,:);

                    for tr = 1:thisDecompAxx.nTrl
                        tR = corr(est_signal(:,:,tr),ref_signal);                       
                        R2 = tR(1:2,:).^2;
                        if R2(1,1)<R2(1,2)% if the order of components is flipped
                            R2 = R2(:,2:-1:1);
                        end
                        trial_R2(:,:,tr) = R2;
                    end             
                    
                    residuals.(this_decomp_method){s}(:,1:2,nLambda_idx,draw_idx) = squeeze(mean(trial_R2,3))'; % average residual over trials
                    clear trial_R2;
                end
                snrs_orig{s}(:,nLambda_idx,draw_idx)=mean(mean(thisAxx.Amp(signal_freq_idxs,:,:).^2)./mean(thisAxx.Amp(noise_freq_idxs,:,:).^2),3);
                  
            end
        end

        %
               
    end
    EEGData = {};
    EEGAxx = {} ;
end

%%
% scalp plots
Subject_idx = 1;
source_pattern = Source_pattern(:,:,Subject_idx);
n_comps = 2;
decomp_methods = {'pca','ssd'};
FigH = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
nCols = 1+n_comps*length(decomp_methods) ;
nRows = length(Lambda_list);
set(FigH,'Unit','Inch','position',[5, 5, 18, 9],'color','w');
rowCounter = 1 ;

for src = 1:2%n_comps
    S = subplot_tight(nRows,nCols,1+nCols*(floor(nRows/2)+src-1),[-0.01,-0.005]);
    set(S,'position',get(S,'position')+[.05 0 0 0])
    mrC.Simulate.PlotScalp(source_pattern(:,src),['Source' num2str(src)]);
    caxis([-max(abs(source_pattern(:,src))) max(abs(source_pattern(:,src)))]);
end

for nLambda_idx = 1:length(Lambda_list)
    colCounter = 1 ;
    for decomp_method_idx=1:length(decomp_methods)
        this_decomp_method = decomp_methods{decomp_method_idx};
        for comp_idx = 1:min(n_comps,    size(A.(this_decomp_method){Subject_idx}{nLambda_idx}{1},2 ))
            Sub = subplot_tight(nRows,nCols,1+(colCounter-1)*n_comps+comp_idx+(rowCounter-1)*nCols,[-0.01,-0.005]);
            set(Sub,'position',get(Sub,'position')+[0.05-(comp_idx*.025) -0.02 0 0]);

            Topo = A.(this_decomp_method){Subject_idx}{nLambda_idx}{1}(:,comp_idx);
            if abs(min(Topo))>max(Topo), Topo = -1*Topo;end
            mrC.Simulate.PlotScalp(Topo);
            if comp_idx ==1
                this_title = sprintf('SNR = %i (dB)',round(10*log10(Lambda_list(nLambda_idx))));
                T = title(this_title);
                set(T,'position',get(T,'position')+[2  0 0]);
            end    

            if rowCounter == 1
    %             axes('position',[.18+(.10*(comp_idx-1))+(.22*(decomp_method_idx-1)) .98 .2 .1 ]);
    %             this_title = sprintf('%s %1i',this_decomp_method,  comp_idx) ;
    %             text(0,0,this_title,'fontsize',12,'fontweight','bold'); axis off;
                if comp_idx == 1
                    axes('position',[.38+(.4*(decomp_method_idx-1))-(.08) .99 .2 .1 ]);
                    this_title = sprintf('%s %i',this_decomp_method,comp_idx) ;
                    text(0,0,this_title,'fontsize',12,'fontweight','bold'); axis off;
                else
                    axes('position',[.38+(.4*(decomp_method_idx-1))+.1 .99 .2 .1 ]);
                    this_title = sprintf('%s %i',this_decomp_method,comp_idx) ;
                    text(0,0,this_title,'fontsize',12,'fontweight','bold'); axis off;
                end
            else
                    this_title = sprintf('%i trials',Lambda_list(nLambda_idx));
            end
            caxis([-max(abs(Topo)) max(abs(Topo))]);
        end
        colCounter = colCounter +1 ;
    end
    rowCounter = rowCounter+1 ;
end
colormap(jmaColors('coolhotcortex'))
set(FigH,'Unit','Inch','position',[5, 5, 10, 7],'color','w');
export_fig(FigH,fullfile(FigPath,'SpatialFilters_topographies_SNR_average'),'-pdf');
%close;

 %%  plots angular error of topographies
snrs_orig = cellfun(@(x) x(1:10,:,:),snrs_orig,'uni',false);
snrs_all = squeeze(mean(cat(4,snrs_orig{:}),3)); 
if false
    snrs_inp = 10*log10(squeeze(mean(mean(snrs_all(:,:,:),1),3)));
else    
    snrs_inp = 10*log10(Lambda_list);
end
 FS = 14;
FIG2 = figure;
subplot(1,3,1)
colors = brewermap(4,'Set2') ;
markers = {'-o',':o'};
for comp_idx = 1:2
for decomp_method_idx=1:2%length(decomp_methods)
    this_decomp_method = decomp_methods{decomp_method_idx};
    err_angles.(this_decomp_method) = cellfun(@(x) x(:,1:2,:,:),err_angles.(this_decomp_method),'uni',false);
    errAng_all = squeeze(mean(cat(5,err_angles.(this_decomp_method){:}),4));
    plot(snrs_inp,squeeze(mean(errAng_all(comp_idx,comp_idx,:,:),4)),markers{comp_idx},'LineWidth',2,'MarkerSize',10,'color',colors(decomp_method_idx,:));

    hold on
end
end


%set(gca,'xtick',10*log10(Lambda_list),'xticklabel',arrayfun(@num2str,round(log10(Lambda_list)*10),'uni',false));
xlim([min(snrs_inp)-.2*(abs(min(snrs_inp))) max(snrs_inp)*1.2]);

xlabel('SNR (dB)')
ylabel('Error Angle')

set(gca,'fontsize',FS)

% plot snrs
subplot(1,3,2)
comp_idx =1;

for comp_idx = 1:2
    for decomp_method_idx=1:2%length(decomp_methods)
    this_decomp_method = decomp_methods{decomp_method_idx};
    snrs.(this_decomp_method) = cellfun(@(x) x(1:2,:,:),snrs.(this_decomp_method),'uni',false);
    snrs_all = squeeze(mean(cat(4,snrs.(this_decomp_method){:}),3));
    plot(snrs_inp,10*log10(squeeze(mean(snrs_all(comp_idx,:,:),3))),markers{comp_idx},'LineWidth',2,'MarkerSize',10,'color',colors(decomp_method_idx,:));
    hold on
    end
end
xlabel('SNR (dB)')
ylabel('Output SNR (dB)')
%set(gca,'xtick',10*log10(Lambda_list),'xticklabel',arrayfun(@num2str,round(log10(Lambda_list)*10),'uni',false));
xlim([min(snrs_inp)-.2*(abs(min(snrs_inp))) max(snrs_inp)*1.2]);
set(gca,'fontsize',FS)

% plot residual
subplot(1,3,3)
for comp_idx = 1:2
    for decomp_method_idx=1:2%length(decomp_methods)
        this_decomp_method = decomp_methods{decomp_method_idx};
        residuals.(this_decomp_method) = cellfun(@(x) x(:,1:2,:,:),residuals.(this_decomp_method),'uni',false);
        residuals_all = squeeze(mean(cat(5,residuals.(this_decomp_method){:}),4));
        plot(snrs_inp,squeeze(mean(residuals_all(comp_idx,comp_idx,:,:),4)),markers{comp_idx},'LineWidth',2,'MarkerSize',10,'color',colors(decomp_method_idx,:));
        hold on
    end
end
%[~, hobj, ~, ~] = legend(decomp_methods(1:2));
[~, hobj, ~, ~] = legend({'pca - comp1','ssd - comp1','pca - comp2','ssd - comp2'},'Location','southeast');
legend('boxoff')
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
ht = findobj(hobj,'type','text');
set(ht,'FontSize',11);
xlabel('SNR (dB)')
ylabel('R2')
set(gca,'fontsize',FS)
%set(gca,'xtick',10*log10(Lambda_list),'xticklabel',arrayfun(@num2str,round(log10(Lambda_list)*10),'uni',false));
xlim([min(snrs_inp)-.2*(abs(min(snrs_inp))) max(snrs_inp)*1.2]);
set(FIG2,'Unit','Inch','position',[5, 5, 18, 5],'color','w');
%export_fig(FIG2,fullfile(FigPath,'SpatialFilters_ErrorPlots_Averaged'),'-pdf');


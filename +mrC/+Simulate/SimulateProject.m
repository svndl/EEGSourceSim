function [EEGData,EEGAxx,EEGData_signal,EEGAxx_signal,sourceDataOrigin,masterList,subIDs,allSubjFwdMatrices,allSubjRois] = SimulateProject(projectPath,varargin)
    
    % Syntax: [EEGData,EEGAxx,sourceDataOrigin,masterList,subIDs] = SimulateProject(projectPath,varargin)
    % Description:	This function gets the path for a mrc project and simulate
    % EEG with activity (seed signal as input) in specific ROIs (input),
    % and pink and alpha noises (noise parameters can be set as input)
    %
    % Syntax:	[EEGData,EEGAxx,sourceDataOrigin,masterList,subIDs] = mrC.RoiDemo(projectPath,varargin)
    % 
%--------------------------------------------------------------------------    
% INPUT:
  % projectPath: a cell string, indicating a  path to mrCurrent project
  % folder with individual subjects in subfolders
    %             
    % 
    %
  %   <options>:
    %
  % (Source Signal Parameters)
    %       signalArray:    a NS x nTrials x seedNum matrix, where NS is the number of
    %                       time samples and seedNum is the number of seed sources
    %                       [NS x 1 x 2 SSVEP sources]
    %
    %       signalsf:       sampling frequency of the input source signal
    %
    %       signalType:     type of simulated signal (visualization might differ for different signals)
    %                       
    %       
    %       signalFF:       a seedNum x 1 vector: determines the fundamental
    %                       frequencis of sources
    %       nTrials:        Number of trials. Noise is redrawn for each trial.
    %
    %       originsource    If true, the original source signal will be
    %                       returned in the output, otherwise it will be
    %                       empty. Be careful about using it, if too many
    %                       trials are requested, you might run out of
    %                       memory.
  
  % (ROI Parameters)
    %       rois            a cell array of roi structure that can be
    %                       extracted from mrC.Simulate.GetRoiList for any
    %                       atlas. The roi structure contain .Type field
    %                       which is the type of atlas used:
    %                       (['func']/'wang'/'glass','kgs','benson')
    %                       the other field is .Name which should be the
    %                       name of ROI and .Hemi which indicates the
    %                       hemisphere, and can be 'L' or 'R' or 'B' for
    %                       both hemisphere.
    %                       Note that the cell array can contain ROIs from
    %                       different atlases
    %                       
    %       roiType:        THIS IS NOT NEEDED IF YOU GIVE THE rois INPUT.  
    %                       string specifying the roitype to use. 
    %                       'main' indicates that the main ROI folder
    %                       /Volumes/svndl/anatomy/SUBJ/standard/meshes/ROIs
    %                       is to be used. (['func']/'wang'/'glass','kgs','benson').
    %
    %
    %       roiSpatfunc     a string indicating which spatial function
    %                       will be used to put the seed signal in ROI
    %                       [uniform]/gaussian
    %       roiSize         number of vertices in each ROI
    %
    %       anatomyPath:  The folder should be for the same subject as
    %                       projectPath points to. It should have ROI forders, default
    %                       cortex file, ..
    
  % (Noise Parameters), all this parameters are defined inside "NoiseParam." structure
    %
    %       mu: This number determines the ratio of pink noise to alpha noise
    %
    %       lambda: This number determines the ratio of signal to noise
    %       
    %       alpha nodes: for now the only option is 'all' which means all visual areas  (maybe later a list of ROIs to put alpha in)
    %
    %       mixing_type_pink_noise: for now only 'coh' is implemented, which is default value
    %
    %       spatial_normalization_type: How to normalize noise and generated signal ['active_nodes']/ 'all_nodes'
    %
    %       distanceType: how to calculate source distances ['Euclidean']/'Geodesic', Geodesic is not implemented yet
    
    

  % (Save Parameters)
    %       Save:           If true, save the simulated data in axx format
    %                       in project folder, for each subject like:
    %                       Projectfolder/nl-00xx/Exp_MATL_HCN_128_Avg/
    %                       The results of all subjects are also saved in a
    %                       file in project folder as Raw_c00x.mat the
    %                       condition number is according to cndNum
    %                       parameter
    %
    %       cndNum:         The condition number for simulated EEG
  
  % (Inverse Parameters) .... should be corrected
    %       inverse:        a string specifying the inverse name to use
    %                       [latest inverse]
    %       doSource:       logical indicating whether to use the inverse to push
    %                       the simulated ROI data back into source space
    %                       true/[false]
    %
% OUTPUT:
    %       EEGData:        a NS x e matrix, containing simulated EEG,
    %                       where NSs is number of time samples and e is the
    %                       number of the electrodes
    %
    %
    %       EEGAxx:         A cell array containing Axx structure of each
    %                       subject's simulated EEG. This output is
    %                       available if the signal type is SSVEP
    %
    %       EEGData_signal: a NS x e matrix, containing simulated EEG
    %                       signal without the noise components,
    %                       where NSs is number of time samples and e is the
    %                       number of the electrodes
    %
    %
    %       EEGAxx_signal:  A cell array containing Axx structure of each
    %                       subject's simulated EEG signal without the
    %                       noise components. This output is
    %                       available if the signal type is SSVEP
    %
    %       sourceDataOrigin: a NS x srcNum matrix, containing simulated
    %                           EEG in source space before converting to
    %                           sensor space EEG, where srcNum is the
    %                           number of source points on the cortical
    %                           meshe
    %
    %       masterList:     a 1 x seedNum cell of strings, indicating ROI names
    %
    %       subIDs:         a 1 x s cell of strings, indicating subjects IDs
    %
%--------------------------------------------------------------------------
 % The function was originally written by Peter Kohler, ...
 % Latest modification: Elham Barzegaran, 03.26.2018
 % Modifications: Sebastian Bosse 8/2/2018
 
%% =====================Prepare input variables============================
 
%--------------------------set up default values---------------------------
opt	= ParseArgs(varargin,...
    'anatomyPath'   , [],...   
    'subSelect'        ,[],...
    'inverse'		, [], ...
    'rois'          , [], ...
    'roiType'       , 'wang',...
    'roiSpatfunc'   , 'uniform',...
    'roiSize'       , 200,...
    'signalArray'   , [],...
    'signalsf'      , 100 ,... 
    'signalType'    , 'SSVEP',...
    'signalFF'      , [],...
    'signalSNRFreqBand' ,[],...
    'NoiseParams'   , struct,...
    'doFwdProjectNoise' ,true, ...
    'OptimizeNoiseParam' ,true,...
    'RedoMixingMatrices', false,...
    'Save'          ,true,...
    'cndNum'        ,1, ...
    'nTrials'       ,1, ...
    'originsource'  ,false...
    );

% Roi Type, the names should be according to folders in (svdnl/anatomy/...)
if ~strcmp(opt.roiType,'main')% THIS SHOUDL BE CORRECTED
    switch(opt.roiType)
        case{'func','functional'} 
            opt.roiType = 'functional';
        case{'wang','wangatlas'}
            opt.roiType = 'wang';
        case{'glass','glasser'}
            opt.roiType = 'glass';
        case{'kgs','kalanit'}
            opt.roiType = 'kgs';
        case{'benson'}
            opt.roiType = 'benson';
        otherwise
            error('unknown ROI type: %s',opt.roiType);
    end
else
end


%------------------set anatomy data path if not defined ---------------------
if isempty(opt.anatomyPath)
    anatDir = getpref('mrCurrent','AnatomyFolder');
    if contains(upper(anatDir),'HEADLESS') || isempty(anatDir) %~isempty(strfind(upper(anatDir),'HEADLESS'))
        anatDir = '/Volumes/svndl/anatomy';
        setpref('mrCurrent','AnatomyFolder',anatDir);
    else
    end
else
    anatDir = opt.anatomyPath;
end

%------------------------Check ROIs class----------------------------------
if ~isempty(opt.rois)
    [opt.rois,FullroiNames,RSubID] = CheckROIsArray(opt.rois);
    if isempty(opt.rois)
        display('Simulation terminated');
        EEGData=[];EEGAxx=[];sourceDataOrigin=[];masterList=[];subIDs=[];
    end
else
    [Roi] = mrC.Simulate.GetRoiClass(projectPathfold,anatDir);
    Roi = cellfun(@(x) x.getAtlasROIs(opt.roiType),Roi,'UniformOutput',false);
end
% -----------------Generate default source signal if not given-------------
% Generate signal of interest
if isempty(opt.signalArray) 
    if isempty(opt.rois)
        [opt.signalArray, opt.signalFF, opt.signalsf]= mrC.Simulate.ModelSeedSignal('signalType',opt.signalType); % default signal (can be compatible with the number of ROIs, can be improved later)
    else 
        [opt.signalArray, opt.signalFF, opt.signalsf]= mrC.Simulate.ModelSeedSignal('signalType',opt.signalType,'signalFreq',round(rand(length(FullroiNames),1)*3+3));
    end
end

if isfield(opt,'signalFF')
    if ~iscolumn(opt.signalFF), opt.signalFF = opt.signalFF';end
end

%% ===========================GENERATE EEG signal==========================
projectPathfold = projectPath;
subIDs = subfolders(projectPathfold,0);
projectPath = subfolders(projectPath,1); % find subjects in the main folder

if isempty(opt.subSelect)
    opt.subSelect = subIDs;
end

Inds = ismember(subIDs,opt.subSelect);
subIDs = subIDs(Inds);
projectPath = cellfun(@(x) fullfile(projectPathfold,x),subIDs,'uni',false);

allFwdMatrices = {} ;
for s = 1:length(projectPath)
    %--------------------------READ FORWARD SOLUTION---------------------------  
    % Read forward
    disp (['Simulating EEG for subject ' subIDs{s}]);
    
    fwdPath = fullfile(projectPath{s},'_MNE_',[subIDs{s}]);
    
    % remove the session number from subjec ID
    SI = strfind(subIDs{s},'ssn');
    if ~isempty(SI)
        subIDs{s} = subIDs{s}(1:SI-2);% -2 because there is a _ before session number
    end
    
    % check if the ROIs and Wang atlas (used for alpha noise) exist for this subject
    alphaRoi = mrC.ROIs([],anatDir);alphaRoi = alphaRoi.loadROIs(subIDs{s},anatDir);
    alphaRoi = alphaRoi.getAtlasROIs('wang');
    
    if sum(strcmp(subIDs{s},RSubID)) || (alphaRoi.ROINum ==0)
        EEGData{s}=[];EEGAxx{s}=[];sourceDataOrigin{s}=[];
        warning(['Skip subject ' subIDs{s} '... ROIs can not be found for this subject! '])
        continue;
    end
     
    if exist([fwdPath '-fwd.mat'],'file') % if the forward matrix have been generated already for this subject
        load([fwdPath '-fwd.mat']);
    else
        fwdStrct = mne_read_forward_solution([fwdPath '-fwd.fif']); % Read forward structure
        % Checks if freesurfer folder path exist
        if ~ispref('freesurfer','SUBJECTS_DIR') || ~exist(getpref('freesurfer','SUBJECTS_DIR'),'dir')
            %temporary set this pref for the example subject
            setpref('freesurfer','SUBJECTS_DIR',fullfile(anatDir,'FREESURFER_SUBS'));% check
        end
        srcStrct = readDefaultSourceSpace(subIDs{s}); % Read source structure from freesurfer
        fwdMatrix = makeForwardMatrixFromMne(fwdStrct ,srcStrct); % Generate Forward matrix
    end
    allSubjFwdMatrices{s} = fwdMatrix ;
    
    % ---------------------------Default ROIs----------------------------------
    seedNum = size(opt.signalArray,2); % Number of seed sources
    
    % Select Random ROIs 
    if isempty(opt.rois)
        % Initialized only for the first subject, then use the same for the rest
        RROI = randperm(Roi{1}.ROINum,seedNum);
        %BE careful about this part %%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        opt.rois = cellfun(@(x) x.selectROIs(RROI),Roi,'UniformOutput',false);% Control for consistency
        [~,M] = max(cellfun(@(x) x.ROINum,opt.rois));
        disp (['Number of ROIs :' num2str(opt.rois{M}.ROINum)]);
        FullroiNames =opt.rois{M}.getFullNames;
        disp(['ROI Names : ' cat(2,FullroiNames{:}) ]);
    end
    masterList = FullroiNames; %cellfun(@(x) [x.Name '_' x.Hemi],opt.rois,'UniformOutput',false); 

%-------------------Generate noise: based on Sebastian's code------------------
    
    % -----Noise default parameters-----
    NS = size(opt.signalArray,1); % Number of time samples
    Noise = opt.NoiseParams;
    Noisefield = fieldnames(Noise);
    % TODO: update tests
    if ~any(strcmp(Noisefield, 'mu')),Noise.mu.pink = 1;Noise.mu.alpha = 1;Noise.mu.sensor = 1;end % power distribution between alpha noise and pink noise ('noise-to-noise ratio')
    if ~any(strcmp(Noisefield, 'lambda')),Noise.lambda = 1/NS/2;end % power distribution between signal and 'total noise' (SNR)
    if ~any(strcmp(Noisefield, 'spatial_normalization_type')),Noise.spatial_normalization_type = 'all_nodes';end% 'active_nodes'/['all_nodes']
    if ~any(strcmp(Noisefield, 'distanceType')),Noise.distanceType = 'Euclidean';end
    if ~any(strcmp(Noisefield, 'Noise.mixing_type_pink_noise')), Noise.mixing_type_pink_noise = 'coh' ;end % coherent mixing of pink noise
    if ~any(strcmp(Noisefield, 'alpha_nodes')), Noise.alpha_nodes = 'all';end % for now I set it to all visual areas, later I can define ROIs for it

    % -----Determine alpha nodes: This is temporary?-----
    %alphaRoiDir = fullfile(anatDir,subIDs{s},'Standard','meshes','wang_ROIs');% alpha noise is always placed in wang ROIs
        alpharoiChunk = alphaRoi.ROI2mat(length(fwdMatrix));
    if strcmp(Noise.alpha_nodes,'all'), Noise.AlphaSrc = find(sum(alpharoiChunk,2)); end % for now: all nodes will show the same alpha power over whole visual cortex  

    disp ('Generating noise signal ...');
    
    % -----Calculate source distance matrix-----
    load(fullfile(anatDir,subIDs{s},'Standard','meshes','defaultCortex.mat'));
    surfData = msh.data; surfData.VertexLR = msh.nVertexLR;
    clear msh;
    
    
    % -----This part calculate mixing matrix for coherent noise-----
    if strcmp(Noise.mixing_type_pink_noise,'coh')
        mixDir = fullfile(anatDir,subIDs{s},'Standard','meshes',['noise_mixing_data_' Noise.distanceType '.mat']);
        
        if ~exist(mixDir,'file') || opt.RedoMixingMatrices% if the mixing data is not calculated already
            spat_dists = mrC.Simulate.CalculateSourceDistance(surfData,Noise.distanceType);
            disp(['Calculating mixing matrix for coherent pink noise ...'])
            noise_mixing_data = mrC.Simulate.GenerateMixingData(spat_dists);
            save(mixDir,'noise_mixing_data','-v7.3');
        else
            disp(['Reading mixing matrix for coherent pink noise ...'])
            load(mixDir,'noise_mixing_data');
        end
    end
    
    band_names = fieldnames(noise_mixing_data.matrices) ;
    nSources = size(noise_mixing_data.matrices.(band_names{1}),2);
    % ----- Generate noise-----
    % calculate coherence in channel space to reduced computational effort
    % in every trial
    for band_idx = 1:length(band_names)
        this_band_name = band_names{band_idx} ;
        C = noise_mixing_data.matrices.(this_band_name); 
        
        % normalize along column (per hemisphere) 
        % see ?Generating nonstationary multisensor signals under a spatial coherence constraint
        %C = C./sqrt(sum((C.^2),1));
        noise_mixing_data.matrices.(this_band_name) = C;
        C_chan = zeros(size(fwdMatrix'));
        for hemi = 1:2 % hemisphere by hemisphere
            source_idxs = (hemi-1)*size(C,1)+1:hemi*size(C,1) ;
            C(:,source_idxs) = C(:,source_idxs)./sqrt(sum((C(:,source_idxs).^2),1));
            C_chan(source_idxs,:) =  C(:,source_idxs)*fwdMatrix(:,source_idxs)'; 
        end
        noise_mixing_data.matrices_chanSpace.(this_band_name) = C_chan;
    end
    
    %noise = zeros(NS, size(fwdMatrix,1), opt.nTrials) ;
    for trial_id =1:opt.nTrials 
        disp(['Trial # ' num2str(trial_id)])
        [PinkNoise(:,:,trial_id),AlphaNoise(:,:,trial_id),SensorNoise(:,:,trial_id)] = mrC.Simulate.GenerateNoise(opt.signalsf, NS, nSources, Noise, noise_mixing_data,Noise.spatial_normalization_type,fwdMatrix,opt.doFwdProjectNoise);   
    end

[noise_sig,Noise,SensorNoise] = mrC.Simulate.FitNoise(opt.signalsf, NS, Noise, PinkNoise,AlphaNoise, SensorNoise,fwdMatrix,opt.doFwdProjectNoise,opt.OptimizeNoiseParam);   

%------------------------ADD THE SIGNAL IN THE ROIs--------------------------
    
    disp('Generating EEG signal ...'); 
 
    subInd = strcmp(cellfun(@(x) x.subID,opt.rois,'UniformOutput',false),subIDs{s});
    allSubjRois{s} = opt.rois{find(subInd)} ;
    [EEGData{s},EEGData_signal{s},sourceData] = mrC.Simulate.SrcSigMtx(opt.rois{find(subInd)},fwdMatrix,surfData,opt,noise_sig,SensorNoise,Noise.lambda,'active_nodes');%Noise.spatial_normalization_type);% ROIsig % NoiseParams
       
    if (opt.nTrials==1) || (opt.originsource) % this is to avoid memory problem
        sourceDataOrigin{s} = sourceData;
    else
        sourceDataOrigin{s} = [];
    end
    
    %visualizeSource(sourceDataOrigin{s}, surfData,opt.signalsf,0)
    %% convert EEG to axx format
    if strcmp(opt.signalType,'SSVEP')
        EEGAxx{s}= mrC.Simulate.CreateAxx(EEGData{s},opt);% Converts the simulated signal to Axx format  
        EEGAxx_signal{s}= mrC.Simulate.CreateAxx(EEGData_signal{s},opt);% Converts the simulated signal to Axx format  
    end
    
%% write output to file 
    
    if (opt.Save)
        SavePath = projectPathfold;
        % prepare mrC simulation project
        if ~exist(fullfile(SavePath,subIDs{s}),'dir')
            mkdir(fullfile(SavePath,subIDs{s}));
        end
        
        % Write axx files
        if ~exist(fullfile(SavePath,subIDs{s},'Exp_MATL_HCN_128_Avg'),'dir')
            mkdir(fullfile(SavePath,subIDs{s},'Exp_MATL_HCN_128_Avg'))
        end
        EEGAxx{s}.writetofile(fullfile(SavePath,subIDs{s},'Exp_MATL_HCN_128_Avg',sprintf('Axx_c0%02d.mat',opt.cndNum)));
        
        % Copy Inverse files
        % copyfile(fullfile(projectPath{s},'Inverses'),fullfile(SavePath,subIDs{s},'Inverses'));
        
        % Write Original source Data
        if ~exist(fullfile(SavePath,subIDs{s},'Orig_Source_Simul'),'dir')
            mkdir(fullfile(SavePath,subIDs{s},'Orig_Source_Simul'));
        end
        SourceDataOrigin = sourceDataOrigin{s};
        save(fullfile(SavePath,subIDs{s},'Orig_Source_Simul',sprintf('Source_c0%02d.mat',opt.cndNum)),'SourceDataOrigin');
    end
end

%% save simulated EEG of all subjects in one file

save(fullfile(projectPathfold,sprintf('SimulatedEEG_c0%02d.mat',opt.cndNum)),'EEGData','EEGAxx','subIDs','masterList');

end

function [ROIsArr,FullroiNames,RSubID] = CheckROIsArray(ROIsArr)
    if sum(abs(diff(cellfun(@(x) x.ROINum,ROIsArr))))~=0
        warning ('Number of ROIs is not the same for all subjects');
        RSub = find(cellfun(@(x) x.ROINum,ROIsArr)==0);
        RSubID = cellfun(@(x) x.subID,ROIsArr(RSub),'UniformOutput',false);
        ROIsArr(RSub)=[];% remove the subjects
        
    else
        RSubID = [];
    end

    [~,M] = max(cellfun(@(x) x.ROINum,ROIsArr));
    disp ('Start simulating EEG...');
    disp (['Number of ROIs :' num2str(ROIsArr{M}.ROINum)]);
    FullroiNames =ROIsArr{M}.getFullNames;
    disp(['ROI Names : ' cat(2,FullroiNames{:}) ]);

    % check the order of ROIs in subjects and make them consistent
    FNames = cellfun(@(x) x.getFullNames,ROIsArr,'UniformOutput',false);
    comps = cellfun(@(x) strcmpi(cat(1,FNames{:}),x), FullroiNames,'UniformOutput',false);
    comps = cat(3,comps{:});comps = sum(comps.*repmat(1:numel(FullroiNames),[size(comps,1) 1 size(comps,3)]),3);
    comps = mat2cell(comps,ones(size(comps,1),1),size(comps,2))';
    ROIsArr = cellfun(@(x,y) x.selectROIs(y),ROIsArr,comps,'UniformOutput',false);
end

function [CrossTalk,Errors,ROISource,ScalpData,LIST,subIDs] = ResolutionMatrices(projectPath,varargin)
    
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

    %
%--------------------------------------------------------------------------
 % Latest modification: Elham Barzegaran, 06.13.2018
 % NOTE: This function is a part of mrC toolboxs

%% =====================Prepare input variables============================
 
%--------------------------set up default values---------------------------
opt	= ParseArgs(varargin,...
    'inverse'		, [], ...
    'subSelect'     ,[],...
    'rois'          , [], ...
    'roiType'       , 'wang',...
    'eccRange'      , [],...
    'figFolder'     , [],...
    'plotting'      , false,...
    'anatomyPath'   , [],...
    'doScalpMap'    ,false,...
    'doScalpNorm'   ,false,...
    'doAUC'         ,false ...
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


%-------Set folder for saving the results if not defined (default is desktop)----------
if isempty(opt.figFolder)
    if ispc,home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    else home = getenv('HOME');end
    opt.figFolder = fullfile(home,'Desktop');
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



%% ===========================GENERATE EEG signal==========================
projectPathfold = projectPath;
projectPath = subfolders(projectPath,1); % find subjects in the main folder
subIDs = subfolders(projectPathfold,0); 
if ~isempty(opt.subSelect)
    Inds = ismember(subIDs,opt.subSelect);
    subIDs = subIDs(Inds);
    projectPath = cellfun(@(x) fullfile(projectPathfold,x),subIDs,'uni',false);
end

if isempty(opt.rois)
    Rois = mrC.Simulate.GetRoiClass(projectPathfold);
else
    Rois = opt.rois;
end

for s = 1:length(projectPath)
    %--------------------------READ FORWARD SOLUTION---------------------------  
    % Read forward
    [~,subIDs{s}] = fileparts(projectPath{s});
    disp (['Simulating EEG for subject ' subIDs{s}]);
    
    fwdPath = fullfile(projectPath{s},'_MNE_',[subIDs{s}]);
    
    % remove the session number from subjec ID
    SI = strfind(subIDs{s},'ssn');
    if ~isempty(SI)
        subIDs{s} = subIDs{s}(1:SI-2);% -2 because there is a _ before session number
    end
    
    % To avoid repeatition for subjects with several sessions
    if s>1
        SUBEXIST = strcmpi(subIDs,subIDs{s});
        if sum(SUBEXIST(1:end-1))==1
            disp('EEG simulation for this subject has been run before');
            continue
        end
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
    

    
    %% Read Inverses
    if ~opt.doScalpMap
        if ~isempty(opt.inverse) 
            invPaths{s} = fullfile(projectPath{s},'Inverses',opt.inverse);
            if exist(invPaths{s},'file')
                curInv = mrC_readEMSEinvFile(invPaths{s});
            else
                error(['Inverse ' opt.inverse ' is not found']);
            end
        else
            warning('Please indicate the inverse name...');
        end
    end
    %

    %% Get the ROIs
    %%%%%%%%%%%%%%%%%%%%% add other atlases later...%%%%%%%%%%%%%%%%%%%%%%%

    subInd = strcmp(cellfun(@(x) x.subID,Rois,'UniformOutput',false),subIDs{s});
    SROI = Rois{find(subInd)};
    if strcmpi(opt.roiType,'benson') && ~isempty(opt.eccRange)
        SROICent = SROI.getAtlasROIs('benson',[0 opt.eccRange(1)]);
        SROISurr = SROI.getAtlasROIs('benson',opt.eccRange);
        SROICS = SROICent.mergROIs(SROISurr);
    else
        SROICS = SROI.getAtlasROIs(opt.roiType);
    end
    [roiChunk, NameList] = SROICS.ROI2mat(size(fwdMatrix,2));
      
    % Make the resolution matrix
    %Resolution = fwdMatrix'*curInv;
    ScalpData{s} = roiChunk.'*fwdMatrix';
    
    if opt.doScalpNorm
        ScalpData{s} = ScalpData{s}./repmat(sum(abs(ScalpData{s}),2),[1 size(ScalpData{s},2)]);
    end
    
    if ~opt.doScalpMap
        % calculate inverse and then normalize
        ROISource{s} = ScalpData{s}*curInv;
        ROISource{s} = ROISource{s}./repmat(max(abs(ROISource{s}),[],2),[1 size(ROISource{s},2)]);
        %------------------------Localization Errors-----------------------
        % These Errors are implemented accroding to Cottereau,B.R., Ales, J.M., Norcia, A.M. Human brain mapping(2012)
        % Relative energy
        E = abs(ROISource{s});for x = 1:size(E,1), E(x,E(x,:)<(max(E(x,:))/5))=0;end% get rid of the noise
        Errors{s}.Relative = sum(E.*roiChunk',2)./sum(E,2);
        
        % Focalization
        Errors{s}.Focalization = sum(((ROISource{s}.*roiChunk')-roiChunk').^2,2)./sum((roiChunk').^2,2);
        
        %AUC
        if opt.doAUC
            % load or calculate distance matrix
%             if ~exist(fullfile(anatDir,subIDs{s},'Standard','meshes' ,'Distance_Euclidean.mat'),'file')
%                 load(fullfile(anatDir,subIDs{s},'Standard','meshes','defaultCortex.mat'));
%                 surfData = msh.data; surfData.VertexLR = msh.nVertexLR;
%                 clear msh;
%                 spat_dists = mrC.Simulate.CalculateSourceDistance(surfData,'Euclidean');
%                 %save(fullfile(anatDir,subIDs{s},'Standard','meshes' ,'Distance_Euclidean.mat'),'spat_dists','-v7.3');
%             else
%                 load(fullfile(anatDir,subIDs{s},'Standard','meshes' ,'Distance_Euclidean.mat'))
%             end
            for r = 1:size(roiChunk,2) % find far and close neighbors and compute AUC
                Roi_verts = find(roiChunk(:,r)>0);
%                 roisize = numel(Roi_verts);
%                 rdists = min(spat_dists(Roi_verts,:));
%                 [~, rdist_ind] = sort(rdists);
%                 
%                 Neibor_c = rdist_ind(roisize+1:2*roisize); % close neibours
%                 Vals = ROISource{s}(r,:);Vals(1:2*roisize)=0;
%                 [Val,Neibor_f] = sort(abs(Vals),'descend');% far neibors
%                 N = 5;
%                 Neibor_f = Neibor_f(1:roisize*N);
%                 
                th = 1:-0.001:0;
                for t = 1:numel(th)
                    TP(r,t) = sum(abs(ROISource{s}(r,Roi_verts))>=th(t));
                    FN(r,t) = sum(abs(ROISource{s}(r,Roi_verts))<th(t));
                    
%                     FP_c = sum(abs(ROISource{s}(r,Neibor_c))>=th(t));
%                     FP_f = sum(abs(ROISource{s}(r,Neibor_f))>=th(t));
                    FP(r,t) = sum(abs(ROISource{s}(r,roiChunk(:,r)==0))>=th(t));
                    TN(r,t) = sum(abs(ROISource{s}(r,roiChunk(:,r)==0))<th(t));
                    
%                     TPR(r,t) = TP(r,t)./roisize;
%                     FPR_c(r,t) = FP_c ./roisize;
%                     FPR_f(r,t) = FP_f ./(roisize*N);
                end
            end
%             Errors{s}.TPR = TPR;
%             Errors{s}.FPR_c = FPR_c;
%             Errors{s}.FPR_f = FPR_f;
            Errors{s}.TP = TP;
            Errors{s}.FN = FN;
            Errors{s}.FP = FP;
            Errors{s}.TN = TN;
              
        end
        
        %---------------------individual figure----------------------------
        INVname = opt.inverse; ind = strfind(INVname,'_');
        if opt.plotting ==1
            figure,
            subplot(2,2,1),mrC.Simulate.VisualizeSourceData(subIDs{s},roiChunk(:,1),anatDir,jmaColors('coolhot')); 
            caxis([-1 1]);title(['V1d 0-2 ' subIDs{s} ' Original']);
            subplot(2,2,2), mrC.Simulate.VisualizeSourceData(subIDs{s},ROISource{s}(1,:),anatDir,jmaColors('coolhot')); 
            Data = ROISource{s}(1,:);
            caxis([-max(Data) max(Data)]);title(['V1d 0-2 ' subIDs{s} ' ' INVname(ind(end)+1:end)]);

            subplot(2,2,3),mrC.Simulate.VisualizeSourceData(subIDs{s},roiChunk(:,13),anatDir,jmaColors('coolhot')); 
            caxis([-1 1]);title(['V1d 2-10 ' subIDs{s} ' Original']);
            subplot(2,2,4), mrC.Simulate.VisualizeSourceData(subIDs{s},ROISource{s}(13,:),anatDir,jmaColors('coolhot')); 
            Data = ROISource{s}(13,:);
            caxis([-max(Data) max(Data)]);title(['V1d 2-10 ' subIDs{s} ' ' INVname(ind(end)+1:end)]);
        end
        %---------------------cross talk matrix----------------------------
        CrossTalk{s} = ROISource{s}*roiChunk./repmat(sum(roiChunk),size(ROISource{s},1),1);
        CrossTalkN{s} = CrossTalk{s}./repmat(max(CrossTalk{s},[],2),[1 length(CrossTalk{s})]);
    else
        CrossTalk{s}=[];
        CrossTalkN{s}=[];
        ROISource{s}=[];
        Errors{s} = [];
    end
    
    LIST = NameList;
end
end


function Inv = mrC_readEMSEinvFile(filename)
% Inv = mrC_readEMSEinvFile(filename)
% returns nChannels x nVertices matrix
% 
% based on emseReadInverse, modified to fread nRows x nCols bytes
% beginning at (assumed) end of header rather than by fseeking back from EOF.
% Thus, this implementation should read .inv files with or without the xml-ish footer.

% $Log: mrC_readEMSEinvFile.m,v $
% Revision 1.3  2009/11/18 01:55:03  nicholas
% merged into main branch
%
% Revision 1.2.2.2  2009/11/12 20:45:18  nicholas
% changed 128 channel check from error to warning
%
% Revision 1.2.2.1  2009/07/21 18:17:31  nicholas
% *** empty log message ***
%
% Revision 1.2  2008/11/07 00:03:00  ales
% Further squashed the emse style inverse binary/text reading bug
%
	fid = fopen(filename,'rb','ieee-le');
	if fid == -1
		error('Error opening %s',filename)
	end
	% Get the magic number
	magicNum = upper(fscanf(fid,'%c',8));
	if strcmp(magicNum,'454D5345') % Magic number is OK, we're reading a real inverse file
		% Next read in the major and minor revs, and other header fields.
		% Based on the file format description in Appendix A of EMSE's help file,
		% we expect exactly ten elements in Header, with the dimensions of the inverse
		% matrix in the 9th and 10th position.  Here's what can go wrong:
		% 1) SSI might revise inverse file header field structure without warning us.
		% 2) There will be two extra fields if "cortical thinning was used", whatever that means.
		% For now, this implementation simply checks whether fscanf returns less than expected number of
		% Header elements, otherwise throwing an error.  It falls to you, dear reader,
		% to implement handling of the remaining possibilities listed above should they ever occur.
		[Header,nHeader] = fscanf(fid,'%d',10);
		% fscanf is not robust to bytes following the last header field that have degenerate ASCII values;
		% so we use fgetl, which seems to behave correctly;
                % These lines did not completely fix the bug. Changed to an explicit fseek,
		% see line below if block: fseek(fid,1,0)
		% nHeader = nHeader + 1;
		% Header(nHeader) = str2num(fgetl(fid));

		if nHeader ~= 10
			error('Expecting 10 header elements, found %d in %s',nHeader,filename)
        	end	
		
		%This line explicity sets the file read position to what we think is the begining of good data.
	        fseek(fid,1,0);
        
		nRows = Header(9);
		nCols = Header(10);
		if nCols ~= 128
% 			error('Expecting 128 columns in inverse, found %d in %s',nCols,filename);
			warning('Expecting 128 columns in inverse, found %d in %s',nCols,filename);
		end
		[Inv,nInv] = fread(fid,nCols*nRows,'float64',0,'ieee-le');
		fclose(fid);
		if nInv ~= nCols*nRows
			error( 'Size of inverse (%d) does not match dimensions in file header (%d*%d=%d)',nInv,nRows,nCols,nRows*nCols)
		end
		Inv = reshape( Inv, nCols, nRows );		% nChannels x nVertices
	else
		error('Magic# in %s = %s, expecting 454D5345',filename,magicNum)
	end
end


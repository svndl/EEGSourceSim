function [fwdMatrices,subIDs] = ReadForwards(projectPath)
% this function reads forward solutions of the subjects in a project

projectPathfold = projectPath;
projectPath = subfolders(projectPath,1); % find subjects in the main folder
subIDs = subfolders(projectPathfold,0); 


%%
for s = 1:numel(projectPath)
    [~,subIDs{s}] = fileparts(projectPath{s});
    
    fwdPath = fullfile(projectPath{s},'_MNE_',[subIDs{s}]);
    
    % remove the session number from subjec ID
    SI = strfind(subIDs{s},'ssn');
    if ~isempty(SI)
        subIDs{s} = subIDs{s}(1:SI-2);% -2 because there is a _ before session number
    end
    
    
    if exist([fwdPath '-fwd.mat'],'file') % if the forward matrix have been generated already for this subject
        load([fwdPath '-fwd.mat']);
        fwdMatrices{s} = fwdMatrix;
    else
        fwdStrct = mne_read_forward_solution([fwdPath '-fwd.fif']); % Read forward structure
        % Checks if freesurfer folder path exist
%         if ~ispref('freesurfer','SUBJECTS_DIR') || ~exist(getpref('freesurfer','SUBJECTS_DIR'),'dir')
%             %temporary set this pref for the example subject
%             setpref('freesurfer','SUBJECTS_DIR',fullfile(anatDir,'FREESURFER_SUBS'));% check
%         end
        srcStrct = readDefaultSourceSpace(subIDs{s}); % Read source structure from freesurfer
        fwdMatrices{s} = makeForwardMatrixFromMne(fwdStrct ,srcStrct); % Generate Forward matrix
    end
end
end
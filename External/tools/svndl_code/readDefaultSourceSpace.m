function [ss] = readDefaultSourceSpace(subjId,sourceId)
%[s] = readDefaultSourceSpace(subjId,)
%
%subjId = string specifying the subject id, e.g. 'skeri0044'
%sourceId = string specifying the source space. Defualt: 'ico-5p'



anatDir = getpref('freesurfer','SUBJECTS_DIR');

%append fs4 if not there
if ~(strncmp(subjId(end-2:end),'fs4',3))
    subjId = [subjId '_fs4'];
end

if ~exist('sourceId','var') || isempty(sourceId),
    sourceId = 'ico-5p';
end


    
ctxFilename=fullfile(anatDir,subjId,'bem',[subjId '-' sourceId '-src.fif']);

if ~exist(ctxFilename,'file')
    error(['Cannot find file: ' ctxFilename]);
end

ss = mne_read_source_spaces(ctxFilename);


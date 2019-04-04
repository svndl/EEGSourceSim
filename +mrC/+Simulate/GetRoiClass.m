function [Roi,subIDs] = GetRoiClass(ProjectPath,anatDir,subSelect)
% Syntax: [Roi,subIDs] = GetRoiClass(projectPath,anatDir)
% Description: Get a project and anatomy path and roi atlas and returns the list of of
%           all ROIs in that atlas, and the list of subjects with this altas ROIs in
%           the project
% INPUT: 
%       ProjectPath:    A link to mrC project
%       anatDir:        The link to the anatomy forlder
%       subSelect (optional): 1 x N array of strings each elements indicates the subject IDs to be
%                            selected
%
% OUTPUT:   
%       Roi:    A NSx1 cell array of ROIs class, where NS is the number of
%       subjects in the project
%       subIDs: A NSx1 array of strings, 
%% Set defaults

if ~exist('anatDir','var') || isempty(anatDir)
    anatDir = getpref('mrCurrent','AnatomyFolder');
end

%% Extract the name of all ROI lists
projectPaths = subfolders(ProjectPath,1); % find subjects in the main folder
subIDs = subfolders(ProjectPath,0);

if exist('subSelect','var') && ~isempty(subSelect)
    Inds = ismember(subIDs,subSelect);
    subIDs = subIDs(Inds);
    projectPaths = cellfun(@(x) fullfile(ProjectPath,x),subIDs,'uni',false);
end


roiList = [];

len = length(projectPaths) ;

for s = 1: len
    
    [~,subIDs{s}] = fileparts(projectPaths{s});
    
    display(['Loading Subject ' num2str(subIDs{s}) ' ROIs']);
    % remove the session number from subjec ID
    SI = strfind(subIDs{s},'ssn');
    if ~isempty(SI)
        subIDs{s} = subIDs{s}(1:SI-2);% -2 because there is a _ before session number
    end
    
    Roi{s} = mrC.ROIs([],anatDir);
    Roi{s} = Roi{s}.loadROIs(subIDs{s},anatDir);
end

end
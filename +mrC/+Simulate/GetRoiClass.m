function [Roi,subIDs] = GetRoiClass(projectPath,anatDir)
% Syntax: [Roi,subIDs] = GetRoiClass(projectPath,anatDir)
% Description: Get a project and anatomy path and roi atlas and returns the list of of
%           all ROIs in that atlas, and the list of subjects with this altas ROIs in
%           the project
% INPUT: 
%       ProjectPath:    A link to mrC project
%       anatDir:        The link to the anatomy forlder
%
% OUTPUT:   
%       Roi:    A NSx1 cell array of ROIs class, where NS is the number of
%       subjects in the project
%       subIDs: A NSx1 array of strings, 
%% Set defaults

if ~exist('anatDir','var') || isempty(anatDir),
    anatDir = getpref('mrCurrent','AnatomyFolder');
end

%% Extract the name of all ROI lists
projectPath = subfolders(projectPath,1); % find subjects in the main folder
roiList = [];

len = length(projectPath) ;

for s = 1: len
    
    [~,subIDs{s}] = fileparts(projectPath{s});
    
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
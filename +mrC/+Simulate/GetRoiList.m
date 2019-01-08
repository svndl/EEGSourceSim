function [Roi,roiList,subIDs] = GetRoiList(projectPath,anatDir,roiType,subject)
% Get a project and anatomy path and roi atlas and returns the list of of
% all ROIs in that atlas, and the list of subjects with this altas ROIs in
% the project

%% Set defaults

if ~exist('anatDir','var') || isempty(anatDir),
    anatDir = getpref('mrCurrent','AnatomyFolder');
end

if ~exist('roiType','var')
    roiType = 'wang';
end

if ~exist('subject','var')
    subject= 'first';
end


%% Extract the name of all ROI lists
projectPath = subfolders(projectPath,1); % find subjects in the main folder
roiList = [];

switch subject
    case 'first'
        len = 1;
    case 'all'
        len = length(projectPath) ;
end

noRoi = zeros(len,1);

for s = 1: len
    [~,subIDs{s}] = fileparts(projectPath{s});
    % remove the session number from subjec ID
    SI = strfind(subIDs{s},'ssn');
    if ~isempty(SI)
        subIDs{s} = subIDs{s}(1:SI-2);% -2 because there is a _ before session number
    end
    
    roiDir = fullfile(anatDir,subIDs{s},'Standard','meshes',[roiType,'_ROIs']);
    if ~exist(roiDir,'dir'),
        error ([roiDir ' dir is not defined for ' subIDs{s}]);
    end
    nroiList =  subfiles(fullfile(roiDir,'/*.mat'),0);
    if nroiList{1}~=0, 
        roiList = unique([roiList; nroiList]); 
    else
        noRoi(s)=1;
    end
end

%% Clean ROI names
Ind1 = strfind(roiList,'_');
Ind2 = strfind(roiList,'-');
roiList_C = (cellfun(@(x,y,z) z(x+1:y-1),Ind1,Ind2,roiList,'UniformOutput',false));
hemi = (cellfun(@(y,z) z(y+1:y+1),Ind2,roiList,'UniformOutput',false));
subIDs=subIDs(noRoi==0);

roiList = unique(cellfun(@(x) x(1:end-4),roiList,'uni',false));

%% Output structure
for r = 1:numel(roiList)
    Roi{r}.Type = roiType;
    Roi{r}.Name = roiList_C{r};
    Roi{r}.Hemi = hemi{r};
end

end
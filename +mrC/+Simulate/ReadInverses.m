function [Inverses,subIDs] = ReadInverses(ProjectPath,InvName,subSelect)

% This function loads in the inverse solutions from a mrC project 

% INPUT:
    % Project Path: The folderwhere the forwards and inverses are stored (mrC project path)
    % InvName: the name of inverse file
    % subSelect (optional): 1 x N array of strings, each elements indicates the subject IDs to be
        % selected
%%   
% Read inverse solutions of a mrC project

    projectPaths = subfolders(ProjectPath,1); % find subjects in the main folder
    subIDs = subfolders(ProjectPath,0);
   
    if exist('subSelect','var') && ~isempty(subSelect)
        Inds = ismember(subIDs,subSelect);
        subIDs = subIDs(Inds);
        projectPaths = cellfun(@(x) fullfile(ProjectPath,x),subIDs,'uni',false);
    end
    
    %%
    if ~exist('InvName','var') || isempty(InvName)
        warning('Indicate the type of inverse')
        Inverses = [];
        return;
    else
        inversename = InvName;
    end
    
    for s = 1:numel(projectPaths)
        invPaths = fullfile(projectPaths{s},'Inverses',inversename);
        if exist(invPaths,'file')
            Inverses{s} = mrC_readEMSEinvFile(invPaths);
        else
            Inverses{s}=[];
        end
    end
end
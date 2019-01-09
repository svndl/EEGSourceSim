function Inverses = ReadInverses(ProjectPath,InvName)
   % Read inverse solutions of a mrC project

    projectPaths = subfolders(ProjectPath,1); % find subjects in the main folder
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

% this function sets up the environment varibales
% please run this before running EEGSourceSim functions

function EEGSourceSimSetUp(varargin)

opt	= ParseArgs(varargin,...
    'EEGSourceSim'   , true,...   
    'Dataset'        ,true,...
    'FieldTrip'		, false ...
    );

    setpref('EEGSourceSim','Version',1.0)
    
    
    %% EEGSourceSim, path
    if opt.EEGSourceSim
        if ~exist('+mrC','dir') || isempty(getpref('EEGSSim','EEGSSimPath'))
            disp('Pick the EEGSourceSim directory');
            EEGSourceSimPath = uigetdir(fileparts(mfilename('fullpath')),'Pick the EEGSourceSim directory');
            if EEGSourceSimPath~=0
                setpref('EEGSSim','EEGSSimPath',EEGSourceSimPath);
                addpath(genpath(EEGSourceSimPath));
            else
                setpref('EEGSSim','EEGSSimPath',[]);
            end
        end
    end
    
    %% Dataset path
    if opt.Dataset
        disp('Pick the Dataset directory');
        DatasetPath = uigetdir(fileparts(mfilename('fullpath')),'Pick the Dataset directory');
        if DatasetPath~=0
            setpref('EEGSSim','DatasetPath',DatasetPath);
            AnatomyPath = fullfile(DatasetPath,'anatomy');
            ProjectPath = fullfile(DatasetPath,'FwdProject');
            if exist(AnatomyPath,'dir') && exist(ProjectPath,'dir')
                setpref('EEGSSim','AnatomyPath',AnatomyPath);
                setpref('EEGSSim','ProjectPath',ProjectPath);
            else
                warning('The ''Dataset'' directory is not correct, please run <a href="matlab:EEGSourceSimSetUp">EEGSourceSimSetUp</a> again')
            end
        else
            setpref('EEGSSim','DatasetPath',[]);
            setpref('EEGSSim','AnatomyPath',[]);
            setpref('EEGSSim','ProjectPath',[]);
        end
    end
    
    %% field Trip path
    if opt.FieldTrip
        disp('Pick the Fieldtrip directory if available, otherwise cancel')
        FTPath = uigetdir(fileparts(mfilename('fullpath')),'Pick the Fieldtrip directory if available, otherwise cancel');
        if (FTPath)~=0
            setpref('EEGSSim','FTPath',FTPath);
            addpath(genpath(FTPath));
        else
            setpref('EEGSSim','FTPath',[]);
        end
    end
end

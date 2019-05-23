
% this script sets up the environment varibales
% please run this before running EEGSourceSim functions

function EEGSourceSimSetUp()

    setpref('EEGSourceSim','Version',1.0)

    % EEGSourceSim, path
    disp('Pick the EEGSourceSim directory');
    EEGSourceSimPath = uigetdir(fileparts(mfilename('fullpath')),'Pick the EEGSourceSim directory');
    if EEGSourceSimPath~=0
        setpref('EEGSSim','EEGSSimPath',EEGSourceSimPath);
        addpath(genpath(EEGSourceSimPath));
    else
        setpref('EEGSSim','EEGSSimPath',[]);
    end

    % fDataset path
    disp('Pick the Dataset directory');
    DatasetPath = uigetdir(fileparts(mfilename('fullpath')),'Pick the Dataset directory');
    if DatasetPath~=0
        setpref('EEGSSim','DatasetPath',DatasetPath);
    else
        setpref('EEGSSim','DatasetPath',[]);
    end

    % field Trip path
    disp('Pick the Fieldtrip directory if available, otherwise cancel')
    FTPath = uigetdir(fileparts(mfilename('fullpath')),'Pick the Fieldtrip directory if available, otherwise cancel');
    if (FTPath)~=0
        setpref('EEGSSim','FTPath',FTPath);
        addpath(genpath(FTPath));
    else
        setpref('EEGSSim','FTPath',[]);
    end
end

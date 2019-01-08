function [WPLI] = ConnectProject(ProjectPath,varargin)

% INPUTS
    % projectPath: Cell array of strings, indicating a list of paths to
    % mrCurrent project folders of individual subjects, these projects
    % should include "AXX TRIALS" matlab exported files

%% Set default values and folders

opt	= ParseArgs(varargin,...
    'plot'		, 0 ...
    );

% find folders of matlab exported Axx files 
AxxPaths = cellfun(@(X) fullfile(X,'Exp_MATL_HCN_128_Avg') ,ProjectPath,'UniformOutput',false); % matlab folders
EXST = cellfun(@(X) exist(X,'dir'), AxxPaths); % Checks if the destination folder exist
AxxPaths = AxxPaths(EXST==7);

load (fullfile(AxxPaths{1},'SsnHeader_ssn'));% Gets number of conditions
nCnd = size(CndParams,1);

%% Read Axx
for s = 1:length(AxxPaths)
    disp(['Reading Subject #' num2str(s) ' data']);
    for cond = 1:nCnd
        AxxStrct = load(fullfile(AxxPaths{s},['Axx_c0' sprintf('%02d',cond) '_trials']));
        EEGAxx{s,cond} = mrC.axx(AxxStrct,0);
        [wpli{s,cond}, csd{s,cond}] = mrC.Connectivity.ConnectAxx(EEGAxx{s,cond});
    end
end

%%

for cond = 1:10
    WPLI(:,:,:,cond) = (mean(cat(4,wpli{:,cond}),4)); % individual connectivity
    PHASM(:,:,:,cond) = (imag(mean(cat(4,csd{:,cond}),4))); % phase difference (average)
    PHASS(:,:,:,cond) = (imag(std(cat(4,csd{:,cond}),[],4)));
    CSD(:,:,:,:,cond) = permute(cat(4,csd{:,cond}),[4 1 2 3]);% group level connevtivity
end

[wplig] = ft_connectivity_wpli(CSD,'dojack',1,'debias',1);
%CSD2 = permute(CSD(:,:,:,:,1:3),[1 5 2 3 4]);
%CSDD = reshape(CSD2,[21 128 128 17]);

%% plots
if opt.plot==1,
    figure,
    for cond = 1:5
        subplot(2,5,cond),imagesc(WPLI(:,:,1,cond));caxis([0 1]);
        subplot(2,5,cond+5),imagesc(WPLI(:,:,1,cond+5));caxis([0 1]);
    end

    figure,
    for cond = 1:5
        subplot(2,5,cond),imagesc(wplig(:,:,1,cond));caxis([0 1]);
        subplot(2,5,cond+5),imagesc(wplig(:,:,1,cond+5));caxis([0 1]);
    end

    M = .9;
    conMap = jmaColors('coolhotcortex');
    figure,
    for cond = 1:5
        subplot(2,5,cond),imagesc((PHASM(:,:,1,cond)));caxis([-M M]);
        subplot(2,5,cond+5),imagesc((PHASM(:,:,1,cond+5)));caxis([-M M]);
    end
    colormap(conMap)


    figure,
    bar(squeeze(PHASM(66,81,1:5,:))');
end

end


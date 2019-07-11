function [wpli, CSDA] = ConnectAxx(EEGAxx,varargin)
% This function get and Axx trial file and calculates connectivity between
% EEG electrodes based on the method that is determined in the input

%% Set default values
opt	= ParseArgs(varargin,...
    'method'        , 'wPLI', ...
    'Permutaion'    , 0,  ...
    'Harmonics'     , 'all' ...
    );
%%
F1 = EEGAxx.i1F1; % indexing in Axx file in according to C language, in matlab 1 should be added

switch opt.Harmonics
    case 'all'
        Freq = F1+1:F1:EEGAxx.nFr;
    case 'even'
        Freq = 2*F1+1:2*F1:EEGAxx.nFr;
    case 'odd'
        Freq = F1+1:2*F1:EEGAxx.nFr;
end

switch opt.method
    case 'wPLI'
        % calculate cross spectrum density matrix with the size nTrl x nCh x nCh x nFr
        EEGSpec = EEGAxx.Cos + 1i*EEGAxx.Sin;
        z = repmat(permute(EEGSpec(Freq,:,:),[3 1 2]),[1 1 1 EEGAxx.nCh]);
        CSD = permute(z .* conj(permute(z,[1 2 4 3])),[1 3 4 2]); 
        [wpli] = ft_connectivity_wpli(CSD,'dojack',1,'debias',1);
        CSDA = squeeze(mean(CSD));
    otherwise
        error('Connectivity method is not recognized...')
end

pval = []; % for now, later permutation will be implemented

end
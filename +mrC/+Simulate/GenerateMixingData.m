function noise_mixing_data = GenerateMixingData(spat_dists, decomp_algo)
% Syntax: noise_mixing_data = GenerateMixingData(spat_dists)
% Description: This function get the spatial distance between sources and loads the
% spatial decay model for coherenc and generates the mixing matrix for noise
% see ﻿10.1121/1.2987429.
% We use the Cholesky decomposition due to complexity considerations. This
% is appropriate as all sources have (approx.) identical PSD and no
% 'perceptual effects' are expected

%% --------------------------------------------------------------------------

if ~exist('decomp_algo','var') || isempty(decomp_algo) 
    decomp_algo = 'cholesky';
end


% preparing mixing matrices
    load('spatial_decay_models_coherence')% this is located in simulate/private folder, it can be obtained by run the code 'spatial_decay_of_coherence.m'

    % calcualting the distances and the coherence takes some time, better to
    % precalculate, write and read
    % mixing_matrices remain constant over all trials!
    band_names = fieldnames(best_model) ;
    noise_mixing_data.band_freqs = band_freqs ;
    noise_mixing_data.mixing_type = 'coh' ;

    % we are assuming isolation between hemispheres: calculate mixing
    % matrices seperately
    
    hWait = waitbar(0,'Calculating mixing matrices ... ');
       
    for freq_band_idx = 1:length(band_names)
        this_band_name = band_names{freq_band_idx} ;
        this_spatial_decay_model = best_model.(this_band_name) ;
        mixing_matrix = zeros(size(spat_dists,2)/2,size(spat_dists,1)) ;
        for hemisphere_idx = 1:2
            this_spat_dists = spat_dists((hemisphere_idx-1)*size(spat_dists,1)/2+(1:size(spat_dists,1)/2),...
                                (hemisphere_idx-1)*size(spat_dists,2)/2+(1:size(spat_dists,2)/2));
            this_coh = this_spatial_decay_model.fun(this_spatial_decay_model.model_params,this_spat_dists);
            this_coh = min(max(this_coh,0),1) ;
                
            if strcmpi(decomp_algo,'cholesky')
                this_mixing_matrix =  chol(this_coh) ;
            elseif strcmpi(decomp_algo,'eigenvalue')
                [V,D] =eig(this_coh);
                this_mixing_matrix = sqrt(D)*V' ;
            else
                error('decomposition method not implemented')
            end
            mixing_matrix(:,(hemisphere_idx-1)*size(spat_dists,1)/2+(1:size(spat_dists,1)/2)) = this_mixing_matrix';
        end
        %noise_mixing_data.band_freqs.(this_band_name) = band_freqs{freq_band_idx} ;
        noise_mixing_data.matrices.(this_band_name) = mixing_matrix;
        waitbar(freq_band_idx/length(band_freqs));
    end       
    close(hWait);
end

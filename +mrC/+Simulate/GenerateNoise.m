function [noise, pink_noise, alpha_noise,sensor_noise] = GenerateNoise(f_sampling, n_samples, n_nodes, NoiseParams, noise_mixing_data, spatial_normalization_type,fwdMatrix)
% Syntax: [noise, pink_noise, pink_noise_uncoh, alpha_noise] = GenerateNoise(f_sampling, n_samples, n_nodes, mu, alpha_nodes, noise_mixing_data, spatial_normalization_type,fwdMatrix)
% Desciption: GENERATE_NOISE Returns noise of unit variance as a combination of alpha
%               activity (bandpass filtered white noise) and spatially coherent pink
%               noise (spectrally shaped white noise)
% INPUT:
    % f_sampling:                   sampling frequency
    % n_samples:                    number of temporal samples to be generated
    % n_nodes:                      number of nodes/vertices to be generated
    % mu:                           power of pink noise/power of alpha activity ('noise-to-noise' ratio)
    % alpha_nodes:                  indices of nodes/vertices carrying alpha activity
    % noise_mixing_data:            data necessary to impose a statistical spatial
    %                               relation on the poink noise
    % spatial_normalization_type:   spatial reference to normalize the
    %                               different noises to
    %                               'all_nodes': normalize to total number of nodes
    %                               'active_nodes': normalize only to nodes where a specific noise has any activity
    % fwdMatrix:                    Matrix of forward model. If given,
    %                               noise is return projected to channel
    %                               space. This reduces computational
    %                               complexity in during noise generation
% OUTPUT:
    % noise: returns a matrix of size [n_samples,n_nodes]
    % pink_noise
    % ..

% Author: Sebastian Bosse
% Latest Modification: EB, 07/17/2018

%% ---------------------------- generate alpha noise------------------------
    %  
    if ~exist('fwdMatrix','var')|| isempty(fwdMatrix)
        doFwdProjection = false;
    else
        doFwdProjection = true ;
    end
    
    alpha_noise = zeros(n_samples,n_nodes);
    alpha_noise(:,NoiseParams.AlphaSrc)  = repmat(GetAlphaActivity(n_samples,f_sampling,[8,12]),[1,length(NoiseParams.AlphaSrc )]); 
    
    if doFwdProjection    
        alpha_noise = alpha_noise*fwdMatrix' ;
    end
    
    if strcmp(spatial_normalization_type,'active_nodes')
        if doFwdProjection
            warning('spatial normalization with active nodes does not really make sense in channel space!')
        end
        n_active_nodes_alpha = sum(sum(abs(alpha_noise))~=0) ;
        alpha_noise = n_active_nodes_alpha*alpha_noise/norm(alpha_noise,'fro') ;
    elseif strcmp(spatial_normalization_type,'all_nodes')
        alpha_noise = alpha_noise/norm(alpha_noise,'fro') ;

    else
        error('%s is not implemented as spatial normalization method', spatial_normalization_type)
    end
    
    
    
%% -----------------------------generate pink noise------------------------
    pink_noise = GetPinkNoise(n_samples, n_nodes);
    % impose coherence on pink noise
    if strcmp(noise_mixing_data.mixing_type,'coh') % just in case we want to add other mixing mechanisms
        % force noise to be spatially coherent within 'hard' frequency
        % ranges
        % for details see: DOI:10.1121/1.2987429
        f = fftshift([-0.5:1/n_samples:0.5-1/n_samples]*f_sampling); % frequncy range
        
        pink_noise_spec = fft(pink_noise,[],1);  
        
        if doFwdProjection
            pink_noise_spec_coh = zeros(size(pink_noise_spec,1), size(fwdMatrix,1)) ;
        else
            pink_noise_spec_coh = zeros(size(pink_noise_spec));
        end
        for band_idx = 1:length(noise_mixing_data.band_freqs)
            % calc coherence for band
            if doFwdProjection
                C = noise_mixing_data.matrices_chanSpace{band_idx}; 
            else
                C = noise_mixing_data.matrices{band_idx}; 
            end
            freq_bin_idxs = (noise_mixing_data.band_freqs{band_idx}(1)<=abs(f))&(abs(f)<noise_mixing_data.band_freqs{band_idx}(2));
            for hemi = 1:2 % hemisphere by hemisphere
                
                if doFwdProjection
                    source_idxs = (hemi-1)*size(C,1)/2+1:hemi*size(C,1)/2 ;
                    pink_noise_spec_coh(freq_bin_idxs,:) = pink_noise_spec_coh(freq_bin_idxs,:)+pink_noise_spec(freq_bin_idxs,source_idxs)*C(source_idxs,:); 
                else
                    source_idxs = (hemi-1)*size(C,2)+1:hemi*size(C,2) ;
                    pink_noise_spec_coh(freq_bin_idxs,source_idxs) =  pink_noise_spec(freq_bin_idxs,source_idxs)*C(source_idxs,:); 
                end
            end
        end
        pink_noise = real(ifft(pink_noise_spec_coh,[],1));
    elseif ~strcmp(noise_mixing_data.mixing_type,'none') 
        error('%s is not implemented as a mixing method',noise_mixing_data.mixing_type)
    end

    if false
       close all
       n_samples = size(pink_noise_spec_coh) ;
       normed_pink_noise_spec_coh = pink_noise_spec_coh/norm(pink_noise_spec_coh,'fro');
       normed_pink_noise_spec = pink_noise_spec*fwdMatrix';
       normed_pink_noise_spec = normed_pink_noise_spec/norm(normed_pink_noise_spec,'fro') ;
       loglog(mean(abs(normed_pink_noise_spec(1:n_samples/2,:)).^2,2),'b')
       hold on
       loglog(mean(abs(normed_pink_noise_spec_coh(1:n_samples/2,:)).^2,2),'r')
       

        
    end
    
    if strcmp(spatial_normalization_type,'active_nodes')
        if doFwdProjection
            warning('spatial normalization with active nodes does not really make sense in channel space!')
        end
        n_active_nodes_pink = sum(sum(abs(pink_noise))~=0) ;
        pink_noise = n_active_nodes_pink * pink_noise/norm(pink_noise,'fro') ;
    elseif strcmp(spatial_normalization_type,'all_nodes')
        pink_noise = pink_noise/norm(pink_noise,'fro') ;
    else
        error('%s is not implemented as spatial normalization method', spatial_normalization_type)
    end
    
    sensor_noise = randn(n_samples, size(fwdMatrix,1)) ;
    sensor_noise = sensor_noise/norm(sensor_noise,'fro'); 
%% --------------------combine different types of noise--------------------
    norm_factor = sqrt(NoiseParams.mu.pink^2+NoiseParams.mu.alpha^2+NoiseParams.mu.sensor^2) ;
    noise = NoiseParams.mu.pink/norm_factor*pink_noise + NoiseParams.mu.alpha/norm_factor*alpha_noise + NoiseParams.mu.sensor/norm_factor*sensor_noise ;
    noise = noise/norm(noise,'fro') ;
%% ---------------------------show resulting noise-------------------------
    if false % just to take a look at the noise components, averaged over all channels for power spectrum
        f = [-0.5:1/n_samples:0.5-1/n_samples]*f_sampling; % frequncy range
        t = [0:n_samples-1]/f_sampling ;
        subplot(4,3,1)
        title('pink noise')
        plot(t, pink_noise(:,1:50:end));xlim([0 2]);
        subplot(4,3,2)
        plot(f, mean(abs(fftshift(fft(pink_noise))),2) );xlim([0 max(f)]);
        subplot(4,3,3)
        loglog(f, mean(abs(fftshift(fft(pink_noise))),2) );xlim([0 max(f)]);
        %ylim([0 .2]);

        subplot(4,3,4)
        title('alpha noise')
        plot(t, alpha_noise(:,1:50:end));xlim([0 2]);
        subplot(4,3,5)
        plot(f, mean(abs(fftshift(fft(alpha_noise))),2) );xlim([0 max(f)]);
        subplot(4,3,6)
        loglog(f, mean(abs(fftshift(fft(alpha_noise))),2) );xlim([0 max(f)]);
        
        %plot(f, abs(fftshift(fft(alpha_noise(:,1:50:end)))));xlim([0 max(f)]);
        %ylim([0 .2]);
        
        
        subplot(4,3,7)
        title('sensor noise')
        plot(t, sensor_noise(:,1:50:end));xlim([0 2]);
        subplot(4,3,8)
        plot(f, abs(fftshift(fft(sensor_noise(:,1:50:end)))));xlim([0 max(f)]);
        subplot(4,3,9)
        loglog(f, mean(abs(fftshift(fft(sensor_noise))),2) );xlim([0 max(f)]);

        
        subplot(4,3,10)
        title(sprintf('combined noise with [%1.2f, %1.2f, %1.2f]', NoiseParams.mu.pink, NoiseParams.mu.alpha, NoiseParams.mu.sensor));
        plot(t, noise(:,1:50:end)); xlim([0 2]);
        subplot(4,3,11)
        plot(f, mean(abs(fftshift(fft(noise))),2) );xlim([0 max(f)]);
        subplot(4,3,12)
        loglog(f, mean(abs(fftshift(fft(noise))),2) );xlim([0 max(f)]);

    end
end

function y = GetAlphaActivity(n_samples,sampling_freq,freq_band)

% generate alpha activity as band-pass filtered white noise
% returns a matrix of size [n_samples,n_nodes]
% Author: Sebastian Bosse 03/2017
%--------------------------------------------------------------------------

if nargin <3
    n_trials = 1 ;
end

% generate white noise
x = randn(n_samples,1);

% bandpass white noise according to alpha band
[b,a] = butter(3, freq_band/sampling_freq*2); 
y = filter(b,a, x); 


%normalize to unit variance
y = y./sqrt(mean(abs(y).^2));
end

function pink_noise = GetPinkNoise(n_samples,n_nodes)

% generate pink noise by shaping white noise
% returns a matrix of size [n_samples,n_nodes]
% Author: Sebastian Bosse 01/2018
%--------------------------------------------------------------------------

    M = n_samples + rem(n_samples,2) ;
    n = 1:M ;
    scalings = sqrt(1./n);
    noise_spec = fft(randn(M,n_nodes)).*scalings' ;
    pink_noise = real(ifft(noise_spec))  ;
    pink_noise = pink_noise(1:n_samples,:) ;
    pink_noise = pink_noise./norm(pink_noise,'fro');
end


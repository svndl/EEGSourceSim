function [noise, NoiseParams,sensor_noise] = FitNoise(f_sampling, n_samples,NoiseParams, pink_noise,alpha_noise, sensor_noise,fwdMatrix,doFwdProjection,optimizeParam)
% Fit the Noise-to-Noise ratio accroding to the real resting EEG data

%% --------------------combine different types of noise--------------------
if optimizeParam
    
    % I do not know why spatial coherent function destroys the shape of 1/f
    % spectrum
    if (n_samples/f_sampling)>1
        load('REC_REO_Averagre_Specs_2sec');
    else
        load('REC_REO_Averagre_Specs_1sec');
    end
    
    AgeRange = [20 60];
    Cond = 'REO';
    
    if strcmp(Cond,'REO')
         PSDEEG = AllSubFFTData(:,2);
    elseif strcmp(Cond,'REC')
         PSDEEG = AllSubFFTData(:,1);
    end
    
    PSDEEG = PSDEEG([SubjectData.Age]>=AgeRange(1) & [SubjectData.Age]<=AgeRange(2));
    ASDEEG = squeeze(sqrt(cat(3,PSDEEG{:})));
    MASDEEG = mean(ASDEEG,2);
    MASDEEG = MASDEEG(2:end);% remove the first frequency, since it is affected by low pass filtering of EEG;

    % 1- calculate sensor space amplitude spectrums of noises
    if ~doFwdProjection
        for tr = 1:size(pink_noise,3)
            Pnoise(:,:,tr) = squeeze(pink_noise(:,:,tr)) * fwdMatrix';
            Anoise(:,:,tr) = squeeze(alpha_noise(:,:,tr)) * fwdMatrix';
        end
    else
        Pnoise = pink_noise;
        Anoise = alpha_noise;
    end
    % same parameters as real EEG data
    if (n_samples/f_sampling)>1
        EPLen = f_sampling * 2; % frequency resolution .5 Hz
    else
        EPLen = f_sampling ; % frequency resolution 1 Hz
    end
    
    % Find the frequency indexes
    NoiseFreq = f_sampling / 2*linspace(0,1,EPLen/2+1);
    FreqI = Freq{1};
    FreqI = FreqI(2:end);% remove the first frequency, since it is affected by low pass filtering of EEG;
    NoiseFreqI = ismember(round(NoiseFreq*10)/10,round(FreqI*10)/10);
    
    % Window def
    win = hann(EPLen); % I add windowing here to be consistent with the real data fft calculation 
    for tr = 1:size(pink_noise,3)
        PNoiseASD(:,:,tr) = abs(fft(Pnoise(:,:,tr).*win));
        ANoiseASD(:,:,tr) = abs(fft(Anoise(:,:,tr).*win));
        SNoiseASD(:,:,tr) = abs(fft(sensor_noise(:,:,tr).*win));
    end
    
    PNoiseASD = mean(mean(PNoiseASD(NoiseFreqI,:,:),2),3);
    ANoiseASD = mean(mean(ANoiseASD(NoiseFreqI,:,:),2),3);
    SNoiseASD = mean(mean(SNoiseASD(NoiseFreqI,:,:),2),3);
    
    % 2- FIND THE BEST PARAMETER FIT
%     prob = optimproblem('ObjectiveSense','min');
%     mu = optimvar('mu',3,1,'LowerBound',0,'UpperBound',5);
    ObjectFunc = @(mu) sqrt(sum((MASDEEG-(mu(1)*PNoiseASD+mu(2)*ANoiseASD+mu(3)*SNoiseASD)).^2));
    muR = fminsearch(@(mu) ObjectFunc(mu),[2 2 2]);
    NoiseParams.mu.pink = muR(1);
    NoiseParams.mu.alpha = muR(2);
    NoiseParams.mu.sensor = muR(3);
    
end

% 3- normalize the noises according to that
norm_factor = sqrt(NoiseParams.mu.pink^2+NoiseParams.mu.alpha^2+NoiseParams.mu.sensor^2) ;
NoiseParams.mu.pink = NoiseParams.mu.pink/norm_factor;
NoiseParams.mu.alpha = NoiseParams.mu.alpha/norm_factor;
NoiseParams.mu.sensor = NoiseParams.mu.sensor/norm_factor;

if doFwdProjection
    noise = NoiseParams.mu.pink/norm_factor*pink_noise + NoiseParams.mu.alpha/norm_factor*alpha_noise + NoiseParams.mu.sensor/norm_factor*sensor_noise ;
else
    noise = NoiseParams.mu.pink/norm_factor*pink_noise + NoiseParams.mu.alpha/norm_factor*alpha_noise ; 
    
    for tr = 1:size(noise,3)% NOT SURE ABOUT THIS
        noisetemp = noise(:,:,tr)*fwdMatrix'+NoiseParams.mu.sensor/norm_factor*sensor_noise(:,:,tr);
        sensor_noise(:,:,tr) = NoiseParams.mu.sensor/norm_factor*sensor_noise(:,:,tr)/norm(noisetemp(:,:,tr),'fro'); 
    end
    
end

for tr = 1:size(noise,3)
    noise(:,:,tr) = noise(:,:,tr)/norm(noise(:,:,tr),'fro') ; % in this case be careful about adding sensor noise later %%%%%%%%%%%%%
end    
    
end



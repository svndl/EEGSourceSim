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
    
%% plot the noises

if false
    plot_noise(pink_noise,PNoiseASD,alpha_noise,ANoiseASD,sensor_noise,SNoiseASD,noise,MASDEEG,muR,FreqI);
end
end

function plot_noise(pink_noise,PNoiseASD,alpha_noise,ANoiseASD,sensor_noise,SNoiseASD,noise,MASDEEG,muR,FreqI)
   FIG = figure;
   LW = 1.5;
   FS = 14;
   el = 29;
   % PINK
   S(1,1) = subplot(4,2,1);
   plot(muR(1)*pink_noise(:,el,1),'k','linewidth',LW-.5);
   ylabel('Amplitude','FontSize',FS)
   set(gca,'xtick',0:300:600,'xticklabel',0:2,'fontsize',FS);
   
   S(1,2) = subplot(4,2,2);
   plot(FreqI,muR(1)*PNoiseASD,'k','linewidth',LW);
   ylabel('ASD','FontSize',FS)
   set(gca,'fontsize',FS);
   xlim([1.5 30])
   ylim([0 2])
   
   % ALPHA
   S(2,1)=subplot(4,2,3);
   plot(muR(2)*alpha_noise(:,el,1),'k','linewidth',LW-.5);
   ylabel('Amplitude','FontSize',FS)
   set(gca,'xtick',0:300:600,'xticklabel',0:2,'fontsize',FS);
   
   S(2,2) = subplot(4,2,4);
   plot(FreqI,muR(2)*ANoiseASD,'k','linewidth',LW);
   ylabel('ASD','FontSize',FS);
   set(gca,'fontsize',FS);
   xlim([1.5 30])
   ylim([0 2])
   
   % Sensor
   S(3,1)=subplot(4,2,5);
   plot(muR(3)*sensor_noise(:,el,1),'k','linewidth',LW-.5);
   ylim([-.1 .1])
   ylabel('Amplitude','FontSize',FS)
   set(gca,'xtick',0:300:600,'xticklabel',0:2,'fontsize',FS);
   
   S(3,2)= subplot(4,2,6);
   plot(FreqI,muR(3)*SNoiseASD,'k','linewidth',LW);
   ylabel('ASD','FontSize',FS)
   set(gca,'fontsize',FS);
   xlim([1.5 30])
   ylim([0 2])
   
   % ALL
   S(4,1) = subplot(4,2,7);
   plot(noise(:,el,1),'k','linewidth',LW-.5);
   xlabel('Time(S)','FontSize',FS);
   ylabel('Amplitude','FontSize',FS);
   set(gca,'xtick',0:300:600,'xticklabel',0:2,'fontsize',FS);
   
   S(4,2) = subplot(4,2,8);
   plot(FreqI,(muR(1)*PNoiseASD+muR(2)*ANoiseASD+muR(3)*SNoiseASD),'k','linewidth',LW); 
   hold on; plot(FreqI,MASDEEG,'-','color','r','linewidth',LW);
   xlim([1.5 30])
   ylim([0 2.5])
   xlabel('Frequency(Hz)','FontSize',FS)
   set(gca,'fontsize',FS);
   ylabel('ASD','FontSize',FS)
   
   L = legend('Simulated Noise','Resting State EEG');
   set(L,'fontsize',FS,'box', 'off');
   
   for r = 1:4
       set(S(r,1),'position',get(S(r,1),'position')+[ -.0 -.0 -.0 -.01]);
       set(S(r,2),'position',get(S(r,2),'position')+[ -.02 -.0 -.0 -.01]);
   end
   
    set(FIG,'paperposition',[1 1 8 8 ]);
    set(FIG,'Unit','Inch','position',[12 5 8 8 ],'color','w');
   
    axes('NextPlot','add','position',[.44 .89 .2 .1]);
    text(.0,.5,'1/f-Activity','fontsize',FS); axis off

    axes('NextPlot','add','position',[.45 .67 .2 .1]);
    text(.0,.5,'\alpha-Activity','fontsize',FS); axis off

    axes('NextPlot','add','position',[.44 .45 .2 .1]);
    h=text(.0,.5,'Sensor Noise','fontsize',FS); axis off

    axes('NextPlot','add','position',[.42 .23 .2 .1]);
    h=text(.0,.5,'Combined Noise','fontsize',FS); axis off
    
    print(fullfile('Figures','Noise_modeling.tif'),'-r300','-dtiff');
    export_fig(FIG,fullfile('Figures','Noise_modeling'),'-pdf');
    close;

end


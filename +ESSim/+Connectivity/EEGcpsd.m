function [CSDmat, COHmat, F] = EEGcpsd(Data,varargin)
%% calculate cross spectrum density for EEG in sensor or ROI space, with pwelch method-> good for broad-band signal

%INPUT:
    % Data: is a [time x ROI(Electrodes) X trails] matrix with the time domain EEG or source space data
    % varargin(optional)
        % SF: sampling frequency
        % Type: Type of spectrum estimation ['fft']/'cpsd'/'shirnkld'
        % winLen: length of Window
        % Nov : number of overlap point for window

        %% parse input and assign default values
opt = ParseArgs(varargin,...
    'SF'        ,[],...
    'Type'   ,'fft',...
    'winLen'    ,[],...
    'Nov'       ,0, ...
    'FreqBand'  ,[0 40] ...
    );
 %% calculate cross spectrum density
 
 for tr = 1:size(Data,3)
     %disp(['Calculating CSD trial #' num2str(tr)]);
     % calculating cross spectrums
     Data1 = squeeze(Data(:,:,tr));
     switch (opt.Type)
         
         case 'cpsd'
             for r1 = 1:size(Data,2)
                 for r2= r1:size(Data,2)
                    [Pxy(:,r1,r2,tr),F] = cpsd(Data1(:,r1),Data1(:,r2),hann(opt.winLen),opt.Nov,opt.winLen,opt.SF);
                 end
             end
             
         case 'fft'
             FEEG = fft(Data1,opt.winLen,1);
             F = opt.SF/2*linspace(0,1,opt.winLen/2+1);
             Pxy(:,:,:,tr) = repmat(conj(FEEG),[1 1 size(FEEG,2)]).* permute(repmat(FEEG,[1 1 size(FEEG,2)]),[1 3 2]);
             
         %case 'shirnkld'
             
         otherwise
             error('This type of spectral estimation is not defined');
     end
     FreqInd = (F>=opt.FreqBand(1) & F<=opt.FreqBand(2));
     % calculating Coherence
     CSDM = Pxy(FreqInd,:,:,tr);
     ASD = arrayfun(@(x) repmat(diag(squeeze(CSDM(x,:,:))),[1 size(CSDM,2)]),1:size(CSDM,1),'uni',false);
     ASD = cat(3,ASD{:});
     Cxy(:,:,:,tr) = CSDM./sqrt(permute(ASD,[3 1 2]).*permute(ASD,[3 2 1]));
end
 %FreqInd = (F>opt.FreqBand(1) & F<opt.FreqBand(2));
 CSDmat = Pxy(FreqInd,:,:,:);
 COHmat = Cxy;
 F = F(FreqInd);
end
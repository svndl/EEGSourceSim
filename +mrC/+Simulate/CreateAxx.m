function EEGAxx= CreateAxx(EEGData,opt)
% syntax: EEGAxx= CreateAxx(EEGData,opt)
% This function converts the simulated EEG into Axx format
% INPUTS:
    % EEGDATA: the output of SrcSigMtx function
    
% OUTPUTS:
    % EEGAxx: a mrC.Axx structure containing EEGData
%==========================================================

% Author: Elham Barzegaran, 3/26/2018
% modified by sebastian bosse 8/4/2018

%%
EEGAxx = mrC.axx();

% Determine the window length for FFT (WLF), and window length in time domain (WLT)
if ~isempty(opt.signalFF)
    if prod(mod(opt.signalFF,1)==0)
        f = opt.signalFF(1);
        for i = 1:numel(opt.signalFF)-1
            f = gcd(f,opt.signalFF(i+1));
        end
        WLT = round((f^-1)*opt.signalsf*2);
    else 
        WLT = (opt.signalFF.^-1)*opt.signalsf*2;%window length for each fundamental frequency: half resolution of fundamentals
        WLT = lcms(round(WLT));% least common multiple
    end
    
    if WLT<(opt.signalsf*2)% find a time window for resolution less than .5 Hz
        WLF = WLT*(WLT\(opt.signalsf*2));
    else 
        WLF = WLT;
    end
else
    WLF = opt.signalsf*2;
    WLT = WLF;
end

% frequency resolution
EEGAxx.dFHz = opt.signalsf/WLF;

% Calculate wave signal
NumEp = floor(size(EEGData,1)/WLT);
EEGAxx.Wave = squeeze(permute(mean(reshape(EEGData(1:NumEp*WLT,:),[WLT NumEp size(EEGData,2) size(EEGData,3)]),2),[1 3 2 4]));

% Time resolution in miliseconds
EEGAxx.dTms = 1000/opt.signalsf;

% calculate spectrum
EEGAxx.nTrl = size(EEGAxx.Wave,3);
NumEp = floor(size(EEGData,1)/WLF);
EEGSPEC = fft(squeeze(permute(mean(reshape(EEGData(1:NumEp*WLF,:),[WLF NumEp size(EEGData,2) size(EEGData,3)]),2),[1 3 2 4])),WLF,1);
f = opt.signalsf*(0:(WLF/2))/WLF; % frequencies
FIND = f<=50; % just keep the spectrum less than 50 Hz
EEGSPEC = EEGSPEC(FIND,:,:);

% EEG spectrum elements
EEGAxx.Cos = real(EEGSPEC);
EEGAxx.Sin = imag(EEGSPEC);
EEGAxx.Amp = abs(EEGSPEC);
EEGAxx.nFr = size(EEGAxx.Cos,1);

% indicate fundamental frequency indexes
if ~isempty(opt.signalFF)
    %[~,~,FFI] = intersect(opt.signalFF,round(f*1000)/1000);
    [~,FFI] = min(abs(repmat(opt.signalFF,[1,length(f)])-repmat(f,[numel(opt.signalFF) 1])),[],2);
    EEGAxx.i1F1 = FFI(1)-1;
    if length(FFI)>1, EEGAxx.i1F2 = FFI(2)-1; end
end
% Other axx parameteres
EEGAxx.DataUnitStr = 'Simulation';
EEGAxx.cndNmb = opt.cndNum;

end

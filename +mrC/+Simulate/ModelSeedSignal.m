function [signalOut, FundFreq, SF] = ModelSeedSignal(varargin)
% This is a trial version, later different source signals will be added (NARMAX,...)
% For now only generates SSVEP like signal

% INPUTS:
    % <options>:
        % signalType: - the model for source signal: [Simple], SSVEP, (NARMAX... for later)
        % sf: - Sampling frequency
        % ns: Number of time samples
        % signalFreq: seedNum x 1 vector - fundamental frequencies of the sources, where srcNum is number of seed sources with different signals
        % signalHarmonic: - cell array of seedNum x 1, each of cell arrays should be a vector indicating amplitude for harmonics in one source 
        % signalPhase: - cell array of seedNum x 1, each of cell arrays should be a vector indicating phase for harmonics in one source 
    
% OUTPUTS    
    % signalOut: - ns x srcNum maxtrix,
    % FundFreq: fundamental frequencies
    % SF: samplng frequency

% Written by Elham Barzegaran, this code uses some part from Sebastian's Bosse codes
% Last modified 3/7/2018

%% Set up default values for signal parameters

opt	= ParseArgs(varargin,...
    'signalType'		, 'SSVEP', ...
    'signalFreq'       , [],...
    'HarmonicAmp'   , [],...
    'HarmonicPhase'      , [],...
    'reliableAmp',       [],...
    'reliablePhase',       [],...
    'nTrials', 1 ,...
    'sf'            , 100,...
    'ns'            , 1000 ...
    );

% (THIS SHOULD BE UPDATED LATER): To allow multiple fundamental freq in one source
if isempty(opt.signalFreq), opt.signalFreq = [2,5];end % initialize fundamental frequency

seedNum = numel(opt.signalFreq);% Number of seed sources

% This part determines the number of harmonics for each source
if strcmp(opt.signalType,'Simple')% just for test, sinusoidal signal
    
   NH = ones(seedNum,1) ; 
   
elseif strcmp(opt.signalType,'SSVEP')% number of harmonics for SSVEP sources, default: 6
    
   if ~(isempty(opt.signalHarmonic)) 
       NH = cellfun('length',opt.signalHarmonic);
   elseif ~(isempty(opt.signalPhase)) 
       NH = cellfun('length',opt.signalPhase); 
   else
    NH = ones(seedNum,1)*6; 
   end
   
end

% If harmonics amplitudes and phases are not defined, they are initialized here randomly!!!
if isempty(opt.signalHarmonic)
    for s = 1:seedNum, opt.signalHarmonic{s} = (rand(1,NH(s))*5); end% random amplitude of harmonics between [0 5]
end
if isempty(opt.signalPhase)
    for s = 1:seedNum, opt.signalPhase{s} = ((rand(1,NH(s))*2-1)*pi); end % random phase of harmonics between [-pi pi]
end  

if isempty(opt.reliableAmp)% all sources are reliable in amplitude
    opt.reliableAmp = ones(1,seedNum);
end

if isempty(opt.reliablePhase)% all sources are reliable in phase
    opt.reliablePhase = ones(1,seedNum);
end

% Checks if the input parameters match
if (seedNum~=numel(opt.signalHarmonic)) || (seedNum~=numel(opt.signalPhase))
    error('Source parameters do not match');
end
for h = 1: numel(opt.signalHarmonic)
    if length(opt.signalHarmonic{h})~=length(opt.signalPhase{h})
        error('Harmonic parameters do not match')
    end
end

%% Generate SSVEP signal
signalOut = zeros(opt.ns,opt.nTrials,seedNum);
t = (0:opt.ns-1)/opt.sf ;


for source_idx = 1:seedNum 
    if opt.reliableAmp(source_idx) && opt.reliablePhase(source_idx) % model reliable sourse
        signalOut(:,:,source_idx) = repmat(sum(opt.signalHarmonic{source_idx}'.*cos(2*pi*opt.signalFreq(source_idx)*[1:length(opt.signalHarmonic{source_idx})]'*t+opt.signalPhase{source_idx}')),[opt.nTrials,1])';
    else % model unreliable sourse
        if opt.reliableAmp(source_idx)
            amps = repmat(opt.signalHarmonic{source_idx},opt.nTrials,1); 
        else
            amps = rand(1,opt.nTrials)' *opt.signalHarmonic{source_idx};
        end
        if opt.reliablePhase(source_idx)
            phs = repmat(opt.signalPhase{source_idx},opt.nTrials,1);
        else
            phs = 2*pi*rand(opt.nTrials,length(opt.signalPhase{source_idx})) ;
        end
        for trial_idx = 1:opt.nTrials
            signalOut(:,trial_idx,source_idx) = sum(amps(trial_idx,:)'.*cos(2*pi*opt.signalFreq(source_idx)*[1:length(opt.signalHarmonic{source_idx})]'*t+phs(trial_idx)')) ;
        end
    end
end
%%
FundFreq = opt.signalFreq;
SF = opt.sf;

%% plot the signal
if false
    figure;
    subplot(1,2,1),plot(1/opt.sf:1/opt.sf:2,signalOut(1:opt.sf*2,:));
    xlabel('Time(s)');
    
    fsig = fft(signalOut,opt.sf);
    f = 0:round(opt.sf/2);
    subplot(1,2,2),bar(f,abs(fsig(1:round(opt.sf/2)+1,:)),2);
    xlabel('Frequency (Hz)'); ylim([0 200]);
end
end


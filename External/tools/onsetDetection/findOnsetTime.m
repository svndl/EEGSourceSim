function [onsetTime,onsetIx] = findOnsetTime(timeSeries,timeUnits,timeLabel,onsetType,varargin)
% [onsetTime,onsetIx] = findOnsetTime(timeSeries,timeUnits,timeLabel,onsetType,['optionalVarName',optionalVarValue,...])
% 
% Returns the onset time and index of an EEG feature in a 1D time series.
%   
% INPUTS:
% timeSeries contains the 1xN vector of data whose onset time is to be 
%   detected 
% timeUnits is a 1xN vector containing the corresponding physical units of
%   time
% timeLabel is a 1xN vector containing -1 where the baseline interval
%   occurs and +1 where the stimulus occurs. The data in timeSeries will be
%   baselined by subtracting the average value during the baseline
%   interval. Baseline interval values are also used by some algorithms to
%   determine thresholds.
% onsetType is a string specifying which algorithm to use:
%
%   'OsmanEtAl1992' estimates noise variability to set a threshold
%   for onset. It requires the absolute signal to remain above threshold 
%   for at least two consequtive 50 msec windows. The additional variable 
%   'allSeries' should be supplied as an optional variable preceded by the 
%   string 'allSeries'. allSeries (MxN) contains M observations that will
%   collectively be used to estimate the noise distribution. Osman et al.
%   (1992) computed the noise estimate using all of the time series going
%   into a condition comparison. If 'allSeries' is not supplied, the noise 
%   estimate will rely on timeSeries. Note, timeUnits are assumed to be in 
%   SECONDS. A different window length can be specified via the optional
%   parameter 'windowLenSec'.
%
%   '2.5stdThresh' is very similar to 'OsmanEtAl1992' except it does not
%   require two following windows to exceed the threshold.
%   
% to suppress plotting, supply the optional parameter 'makePlot' with value
% false or 0.
%
% OUTPUTS:
% onsetTime in the units of timeUnits
% onsetIx is the index in timeSeries corresponding to onsetTime


if nargin<4
    onsetType = 'OsmanEtAl1992';
end

makePlot = true;
printInfo = true;
signedOnset = false;
windowLenSec = nan;
hardThresh = nan;
areaFraction = nan;
windowLogical = nan;
allSeries = nan;
peakFraction = nan;

if nargin>4 % optional parameters were specified
    for k = 1:2:length(varargin)
        switch varargin{k}
            case 'allSeries'
                allSeries = varargin{k+1};
            case 'hardThresh'
                hardThresh = varargin{k+1};
            case 'windowLenSec'
                windowLenSec = varargin{k+1};
            case 'makePlot'
                makePlot = varargin{k+1};
            case 'printInfo'
                printInfo = varargin{k+1};
            case 'areaFraction'
                areaFraction = varargin{k+1};
            case 'windowLogical'
                windowLogical = logical(varargin{k+1});
            case 'signedOnset'
                signedOnset = varargin{k+1};
        end
    end
end

onsetTime = nan;
onsetIx = nan;

% data must always be baselined beforehand
baselineInterval = timeLabel == -1;
stimInterval = timeLabel == 1;
timeSeries = baselineData(timeSeries, baselineInterval);
if ~isnan(allSeries)
    allSeries = baselineData(allSeries, baselineInterval);
end

switch onsetType
    case 'OsmanEtAl1992'
        if printInfo, fprintf('Using the method described by Osman et al. (1992) ...\n'); end
        if isnan(allSeries)
            % use timeSeries baseline interval as noise distribution
            if printInfo, fprintf('Warning: no distribution provided for noise estimate, using timeSeries.\n'); end
            allSeries = timeSeries;          
        end
        if isnan(windowLenSec)
            windowLenSec = 0.05;
        else
            if printInfo, fprintf('Warning: you are using a window length of %.0f msec instead of the standard 50 msec.\n',windowLenSec.*1000); end
        end
        windowLen = find(abs(timeUnits - (timeUnits(1)+windowLenSec))<1e-2, 1 );
        noiseDist = allSeries(:,baselineInterval);
        noiseStd = sqrt(nanvar(noiseDist(:)));
        criterion = noiseStd * 2.5;
        stimSeries = timeSeries(stimInterval);
        timeStim = timeUnits(stimInterval);
        if signedOnset
            if signedOnset == 1
                ix = find(stimSeries>criterion);
            elseif signedOnset == -1
                ix = find(stimSeries<criterion*signedOnset);
            else
                error('signedOnset must be either 1 or -1');
            end
        else
            ix = find(abs(stimSeries)>criterion);
        end
        overThresh = false;
        k = 1;
        while k<length(ix) && ~overThresh && ix(k)<(length(stimSeries)-2*windowLen)
            win1ix = ix(k):ix(k)+windowLen-1;
            win1mean = nanmean(abs(stimSeries(win1ix)));
            win2ix = ix(k)+windowLen:ix(k)+2*windowLen-1;
            win2mean = nanmean(abs(stimSeries(win2ix)));
            if win1mean>criterion && win2mean>criterion
                overThresh = true;
            else 
                k = k + 1;
            end
        end
        if overThresh
            onsetTime = timeStim(ix(k));
            onsetIx = ix(k)+sum(baselineInterval);
        end
    case '2.5stdThresh'
        if printInfo, fprintf('Using a single threshold based on 2.5*noise standard deviation ...\n'); end
        if isnan(allSeries)
            % use timeSeries baseline interval as noise distribution
            if printInfo, fprintf('Warning: no distribution provided for noise estimate, using timeSeries.\n'); end
            allSeries = timeSeries;         
        end
        noiseDist = allSeries(:,baselineInterval);
        noiseStd = sqrt(nanvar(noiseDist(:)));
        criterion = noiseStd * 2.5;
        stimSeries = timeSeries(stimInterval);
        timeStim = timeUnits(stimInterval);
        if signedOnset
            if signedOnset == 1
                ix = find(stimSeries>criterion,1);
            elseif signedOnset == -1
                ix = find(stimSeries<criterion*signedOnset,1);
            else
                error('signedOnset must be either 1 or -1');
            end
        else
            ix = find(abs(stimSeries)>criterion,1);
        end
        if ~isempty(ix)
            onsetTime = timeStim(ix);
            onsetIx = ix+sum(baselineInterval);
        end
    case 'HardThresh'
        if isnan(hardThresh)
            hardThresh = 1;
        end
        if printInfo, fprintf('Using a single hard threshold = %1.2f ...\n',hardThresh); end
        criterion = hardThresh;
        stimSeries = timeSeries(stimInterval);
        timeStim = timeUnits(stimInterval);
        ix = find(abs(stimSeries)>criterion,1);
        if ~isempty(ix)
            onsetTime = timeStim(ix);
            onsetIx = ix+sum(baselineInterval);
        end
    case 'FractionalArea'
        if isnan(areaFraction)
            areaFraction = 0.25;           
        end
        if printInfo, fprintf('Finding point at which %2.2f%% of the area is reached ...\n',areaFraction.*100); end
        if isnan(windowLogical)
            windowLogical = stimInterval;
            if printInfo, fprintf('Using default window of entire stimulus interval from %1.4f to %1.4f time units.\n',timeUnits(find(windowLogical==1,1)),timeUnits(find(windowLogical==1,1,'last'))); end
        else
            if printInfo, fprintf('Using supplied window from %1.4f to %1.4f time units.\n',timeUnits(find(windowLogical==1,1)),timeUnits(find(windowLogical==1,1,'last'))); end
        end
        dataToSum = timeSeries(windowLogical);
        timeStim = timeUnits(windowLogical);
        dataCum = cumsum(abs(dataToSum));
        dataCum = dataCum./sum(abs(dataToSum));
        ix = find((dataCum-areaFraction)>=0,1);
        onsetTime = timeStim(ix);
        onsetIx = ix+find(windowLogical==1,1)-1;     
    case 'FractionalPeakLatency'
        if isnan(peakFraction)
            peakFraction = 0.5;
            if printInfo, fprintf('Finding 1st point at which at least %2.2f%% of the absolute maximum value is reached ...\n',peakFraction.*100); end
        end
        if isnan(windowLogical)
            windowLogical = stimInterval;
            if printInfo, fprintf('Using default window of entire stimulus interval from %1.4f to %1.4f time units.\n',timeUnits(find(windowLogical==1,1)),timeUnits(find(windowLogical==1,1,'last'))); end
        else
            if printInfo, fprintf('Using supplied window from %1.4f to %1.4f time units.\n',timeUnits(find(windowLogical==1,1)),timeUnits(find(windowLogical==1,1,'last'))); end
        end
        if printInfo, fprintf('NOTE! This method does not use interpolation, so onset times will be overestimated for low resolution data.\n'); end
        if printInfo, fprintf('NOTE! It is suggested to use high resolution data when possible.\n'); end
        dataToScan = timeSeries(windowLogical);
        timeStim = timeUnits(windowLogical);
        [peakVal,peakIx] = max(abs(dataToScan));
        criterion = peakVal*peakFraction;
        ix = find(abs(dataToScan(1:peakIx))>=criterion,1);
        onsetTime = timeStim(ix);
        onsetIx = ix+find(windowLogical==1,1)-1;
        if printInfo, fprintf('Detected max of %2.2f (criterion = %2.2f) at %2.4f time units\n',peakVal,criterion,timeStim(peakIx)); end
        if printInfo, fprintf('Detected onset at %2.4f time units (value = %2.2f).\n',onsetTime,dataToScan(ix)); end
    case 'SlopeToPeakZeroInterp'
        % how?
    case 'DoubleThresh'
        % thresh & derivative thresh        
end

if makePlot
    figure;
    set(gca,'FontSize',14)
    hold on;
    plot(timeUnits,zeros(1,length(timeUnits)),'k-')
    plot([0 0],[floor(min(timeSeries)) ceil(max(timeSeries))],'k-')
    if strcmp(onsetType,'HardThresh')
        plot(timeUnits,hardThresh.*ones(1,length(timeUnits)),'r--')
        plot(timeUnits,-hardThresh.*ones(1,length(timeUnits)),'r--')
    end
    plot(timeUnits,timeSeries,'b-','LineWidth',1.5)
    plot(timeUnits(baselineInterval),timeSeries(baselineInterval),'k','LineWidth',1.5)
    grid on;
    if ~isnan(onsetTime)
        plot([onsetTime onsetTime],[floor(min(timeSeries)) ceil(max(timeSeries))],'r--','LineWidth',1.5)
        text(onsetTime+0.025,min(timeSeries),sprintf('t = %2.2f',onsetTime),'Color','r','FontSize',14)
    else
        text(min(timeUnits)+0.05,floor(min(timeSeries))+0.5,'Onset not detected.','Color','r','FontSize',14)
    end
    if strcmp(onsetType,'FractionalArea')
        xx = timeUnits(windowLogical);
        yy = timeSeries(windowLogical);
        h = fill([xx xx(end:-1:1)],[yy zeros(1,length(xx))],'r');
        set(h,'LineStyle','none','FaceAlpha',0.2);
    end
    title(['Onset detection using method: ', onsetType])
    xlabel('Time')
    ylabel('Baselined Data Value')
end

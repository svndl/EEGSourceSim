function [bslnd,meanAmp] = baselineData(data,baseIx)
% [bslnd] = baselineData(data,baseIx)
%
% remove the average value of the datapoints in baseIx from the data and
% return in bslnd. Each row of data will be baselined separately.
%
% data should be a time series supplied with separate observations as
% separate rows.

meanAmp = nan(1,size(data,1));

for s = 1:size(data,1)
    meanAmp(s) = nanmean(data(s,baseIx));
    bslnd(s,:) = data(s,:) - meanAmp(s);
end
    
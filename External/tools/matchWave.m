function [scaleVal,offVal,dataOut] = matchWave(dataA,dataB,doScale,doOffset,matchSep)
    if nargin < 2
        error('function requires to waveforms as input arguments');
    else
    end
    if nargin < 3
        doScale = true;
    else
    end
    if nargin < 4
        doOffset = true;
    else
    end
    if nargin < 5
        matchSep = false;
    else
    end
    
    if ~matchSep % do not do subject by subject matching
        w1 = nanmean(dataA,2);
        w2 = nanmean(dataB,2);
    else
        w1 = dataA;
        w2 = dataB;
    end
    for n = 1:size(w2,2)
        scaleVal(n) = 1;
        offVal(n) = 0;
        % do scaling first, then offset
        if doScale
            scaleVal(n) = fminsearch(@(x) sum ( abs (w1(:,n)-w2(:,n).*x)),1);
        else
        end
        if doOffset
            offVal(n) = fminsearch(@(x) sum ( abs (w1(:,n)-w2(:,n).*scaleVal(n)+x)),0);
        else
        end
        if matchSep
            dataOut(:,n) = dataB(:,n).*scaleVal(n)+offVal(n);
        else
            dataOut = dataB.*scaleVal(1)+offVal(1); % only compute output data once
        end
    end
end

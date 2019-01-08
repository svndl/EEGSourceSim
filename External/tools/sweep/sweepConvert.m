function [sweepValue,respTimes] = sweepConvert( coreDuration,sweepBounds,respTimes,sweepType )
    if nargin < 4 || isempty(sweepType)
        sweepType = 'log';
    else
    end
    if nargin < 3 || isempty(respTimes)
        respTimes = linspace(0,coreDuration,101);
        respTimes = respTimes(2:end);
    else
    end
      
    % first validate the sweep bounds
    if strcmp( sweepType,'log' )
        sweepSign = sign(sweepBounds(1)); % obtain sign for use with negative log sweeps.
        if sweepSign == 0
            error('sweep start value cannot be == 0');
        else
        end
        sweepBounds = abs(sweepBounds);
        if sweepBounds(1) < 0.0001
            warning('sweep start value lower than 0.0001, set to 0.0001');
        else
        end
        if sweepBounds(2) < 0.0001
            warning('sweep end value lower than 0.0001, set to 0.0001');
        else
        end
        sweepBounds = log(sweepBounds);
    else
    end
    
    % compute the sweep fraction
    sweepFraction = respTimes./coreDuration;
    
    % and the values
    sweepValue = sweepBounds(1) + sweepFraction * ( sweepBounds(2) - sweepBounds(1) );
    if strcmp( sweepType , 'log' )
        % convert back to non-log value
        sweepValue = sweepSign * exp( sweepValue); % restore sign for negative log sweep
    else
    end
end
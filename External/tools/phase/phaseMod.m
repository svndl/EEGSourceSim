function outPhase = phaseMod(inPhase,baseShift,convertToDeg)
    if nargin < 3
        convertToDeg = false;
    else
    end
    if nargin < 2
        baseShift = pi/2;
    else
        if baseShift > 2*pi
            error('baseShift must be given in radians');
        else
        end
    end
    outPhase = 2*pi-mod(inPhase+baseShift,2*pi);
    
    %if outPhase < 0 
    %    disp('dude');
    %    outPhase(outPhase < 0) = 2*pi-outPhase(outPhase < 0);
    %else
    %end
    
    if convertToDeg
        outPhase = outPhase.*180/pi;
    else
    end
    
end


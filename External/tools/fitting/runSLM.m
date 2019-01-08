function [slmOut,slmOnset,slmOffset,slmSlope] = runSLM(xdata,ydata,prescription)
    slmOut = slmengine(xdata,ydata,prescription);
    slmOnset = slmOut.knots(2);
    if size(slmOut.knots,1) > 2
        slmOffset = slmOut.knots(3);
        slmSlope = (slmOut.coef(2)-slmOut.coef(3))/(slmOut.knots(2)-slmOut.knots(3));
    else
    end
end
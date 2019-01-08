function [permP,meanAmpRand,meanAmpOrig] = vecPermuteTest(xyData,varargin)
    % Syntax: [permP,meanAmpRand,meanAmpOrig] = vecPermuteTest(xyData,varargin)
    %
    % Returns the results of doing permutation testing of vector values
    %
    % In:
    %   xyData - n x 2 matrix containing n data sample
    %            Function assumes that the 2D data in xyData(:,:,1)
    %            and the optional xyData(:,:,2) are Fourier 
    %            coefficients organized with the real coefficients 
    %            in the 1st column and the imaginary
    %            coefficients in the 2nd column. Rows = samples.
    % 
    % <options>:
    %   nIter - 1-element vector specifying the number of permutations (10,000) 
    %
    % Out:
    %
    %   permP - permutation-based p-value

    nanIdx = sum(isnan(xyData),2)>0;
    xyData = xyData(~nanIdx,:);
    nSubs = size(xyData,1);
    nIter = 10000;
    origAmp = sqrt(xyData(:,1).^2 + xyData(:,2));
    
    % assign phase randomly between 0 and 2pi (360)
    randPhase = rand(nSubs,nIter).*2.*pi;
    randComplex = repmat(origAmp,1,nIter).*exp(1i*randPhase);
    % average over subjects, and convert to amplitude
    realRand = mean(real(randComplex),1);
    imagRand = mean(imag(randComplex),1);
    meanAmpRand = sqrt(realRand.^2+imagRand.^2);
    % get amplitude from original data
    meanAmpOrig = sqrt(mean(xyData(:,1),1).^2+mean(xyData(:,2),1).^2);
    permP = ( nIter - length(find(meanAmpOrig>meanAmpRand)) )/nIter;
end
    
    
    
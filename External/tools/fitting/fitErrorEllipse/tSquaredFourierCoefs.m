function [results] = tSquaredFourierCoefs(xyData,varargin)
    % Syntax: [results] = tSquaredFourierCoefs(xyData,testMu,alphaVal)
    %
    % Returns the results of running Hotelling's t-squared test that the mean
    % of the 2D data in xyData is the same as the mean specified in testMu at
    % the specified alphaVal (0-1).
    %
    % Based on Anderson (1984) An Introduction to Multivariate Statistical 
    % Analysis, Wiley
    %
    % In:
    %   xyData - n x 2 x q matrix containing n data samples
    %            if q = 1, data will be tested against zero    
    %            if q = 2, a paired t-test will be performed,
    %            testing xyData(:,:,1) against xyData(:,:,2). 
    %            2 is the maximum length of the third dimension. 
    %            Function assumes that the 2D data in xyData(:,:,1)
    %            and the optional xyData(:,:,2) are Fourier 
    %            coefficients organized with the real coefficients 
    %            in the 1st column and the imaginary
    %            coefficients in the 2nd column. Rows = samples.
    % 
    % <options>:
    %   testMu - 2-element vector specifying the mean to test against ([0,0]) 
    % 
    %   alphaVal - scalar indicating the alpha value for testing ([0.05])
    %
    %   pdStyle - do PowerDiva style testing (true/[false])
    %
    % Out:
    %
    %   results - struct containing the following fields:
    %             alpha, tSqrdCritical, tSqrd, pVal, H (0 if can't reject null; 1 if
    %             rejected the null hypothesis)

    opt	= ParseArgs(varargin,...
        'testMu'		, [0,0], ...
        'alphaVal'		, 0.05	, ...
        'pdStyle',        false ...
        );

    dims = size(xyData);
    if dims(1) < 2
        error('input data must contain at least 2 samples');
    end
    if dims(2) ~= 2
        error('input data must be a matrix of 2D row samples');
    end
    twoSamples = true;
    if length(dims) < 3 % if no third dimension
        xyData(:,:,2) = zeros(size(xyData));
        twoSamples = false;
    elseif dims(3) > 2
         error('length of third dimension of input data may not exceed two')
    else
    end

    if length(opt.testMu) ~= 2
        error('testMu should be a 2D vector');
    end
    
    % between condition differences
    xy_diff = xyData(:,:,1) - xyData(:,:,2);
    % remove NaNs
    notNaNs = ~any(isnan(xy_diff),2);
    xy_diff = xy_diff(notNaNs,:);
    xyData = xyData(notNaNs,:,:);
    
    N = size(xy_diff,1);
    try
        [sampMu,sampCovMat] = eigFourierCoefs(xy_diff);
    catch
        fprintf('The covariance matrix of xyData could not be calculated, probably your data do not contain >1 sample.');
    end
        
    if opt.pdStyle
        sampCovMat = eye(2);
    else
    end

    results.alpha = opt.alphaVal;

    p = 2; % number of variables
    df1 = p;  %Numerator degrees of freedom.
    df2 = N-p;  %Denominator degrees of freedom.
    
    % Eqn. 2 in Sec. 5.3 of Anderson (1984):
    t0Sqrd = ( (N-1) * p )/ ( N-p ) * finv( 1-opt.alphaVal, p, N - p ); 
    results.tSqrdCritical = t0Sqrd;
    
    try
        invSampCovMat = inv(sampCovMat);    
        % Eqn. 2 of Sec. 5.1 of Anderson (1984):
        tSqrd = N * (sampMu - opt.testMu) * invSampCovMat * (sampMu - opt.testMu)'; 

        tSqrdF = (N-p)/(p*(N-1)) * tSqrd; % F approximation 
        pVal = 1 - fcdf(tSqrdF, df1, df2);  % compute p-value
    catch
        fprintf('inverse of the sample covariance matrix could not be computed.')
        return
    end
    if twoSamples
        % standard deviations for real and imaginary values
        xy_std = std(xyData,0,1);
        % univariate effect size for real and imaginary
        % standardize by average standard deviation over condtions
        % (d_av, equation 10, Läkens, Frontiers in Psych, 2013)
        xy_d = sampMu./sum(xy_std,3);
        
        % compute covariance between real and imaginary values, 
        % and standardize by standard deviation
        % then average over conditions
        % (Del Giudice, personal communication, 2018)
        corMat = (cov([xyData(:,1,1),xyData(:,2,1)])/(xy_std(:,1,1)*xy_std(:,2,1)) + ...
                  cov([xyData(:,1,2),xyData(:,2,2)])/(xy_std(:,1,2)*xy_std(:,2,2)))/2;

        corMat(eye(2)==1) = 1;
        % Mahalanobis' D (equation 1, Del Giudice, Evo. Psych., 2009)
        mahaD = sqrt( xy_d * inv(corMat) * xy_d' );
        % overlapping coefficient (equation 2, Del Giudice, Evo. Psych., 2009)
        OVL = 2*normcdf(-mahaD/2);
        % Cohen?s index of nonoverlap (equation 3, Del Giudice, Evo. Psych., 2009)
        nonOverlap = 1-OVL/(2-OVL); 
    else
        % note that overlap indices do not apply here
        mahaD = sqrt(tSqrd/N);
    end
    results.tSqrd = tSqrd;
    results.pVal = pVal;
    results.H = tSqrd >= results.tSqrdCritical;
    results.df1 = df1;
    results.df2 = df2;
    results.testAmp = abs(complex(sampMu(1),sampMu(2)));
    results.mahalanobisD = mahaD;
    if twoSamples
        results.cohenNonOverlap = nonOverlap;
end


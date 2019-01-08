function [realT,realP,corrT,critVal,clustDistrib]= ttest_permute( inData,maxPerms,timeMask,permMatrix,deletePool,tParams )
        % [realT,realP,corrT,critVal,clustDistrib]= ttest_permute( inData,maxPerms,timeMask,permMatrix,deletePool,tParams )
        % 
        % input:
        %   inData: the difference waveform to test against zero, as a t x s matrix, where t is time and s is subjects.
        % 
        % optional:
        %   maxPerms (scalar): number of permutations to run, if not given every possible permutation will be run
        %   timeMask (1xt logical): true indicates timepoints that should be included in the test, 
        %                           default is to include everything.
        %   permMatrix: matrix of permutations to run, 
        %               useful if multiple tests will be run on the same set of subjects.
        %   deletePool (logical true/[false]): if true, the parallel pool object used to run the tests will be deleted afterwards. Useful if you know you are running the last instance of the function.
        %   tParams (cell): t-test parameters, default is {'dim',2,'alpha',0.05}.
        % 
        % output:
        %   realT (1xt logical): significance of ordinary t-test
        %   realP (1xt double): p-values of ordinary t-test
        %   corrt (1xt logical): significance of corrected t-test
        %   critVal (scalar): critical number of significant time-points
        %                     required for surviving correction
        %   clustDistrib: distribution of significance run lengths for
        %                 permuted data.
        
        % take out NaN data and count subjects
        notNaN = sum(isnan(inData),1)==0;
        inData = inData(:,notNaN);
        numSubs = size(inData,2);
        
        numPerms = 2 ^ numSubs;
        if nargin < 6 || isempty(tParams)
            tParams = {'dim',2,'alpha',0.05};
        else
        end
        if nargin < 5 || isempty(deletePool)
            deletePool = false;
        else
        end
        if nargin < 4 || isempty(permMatrix)
            permMatrix = [];
        end
        if nargin < 3 || isempty(timeMask)
            timeMask = 1:size(inData,1);
        else
        end
        if nargin < 2 || isempty(maxPerms)
            maxPerms = numPerms;
            iPerms = 1:maxPerms;
        else
            % choose a random set of the possible perms/combs
            rand('state',sum(100*clock));
            iPerms = round( numPerms * rand( maxPerms, 1 ) );
        end

        
        if isempty(permMatrix)
            permMatrix = num2bits( iPerms , numSubs );
        else
        end
        
        realData = inData(timeMask,:);
        realT = zeros(size(inData,1),1);
        realP = ones(size(inData,1),1);
        [realT(timeMask),realP(timeMask),CI,stats] = ttest(realData,0,tParams{:});
        clustDistrib = zeros(maxPerms,1);
        poolobj = gcp;
        parfor p = 1:maxPerms % all possible permutations
            fakeData = realData;
            fakeData = fakeData.*repmat(permMatrix(p,:),size(fakeData,1),1);
            %for z=1:size(fakeData,2); fakeData(:,z)=fakeData(:,z).*allPerms(:,p);end;
            [fakeT,P,CI,stats] =  ttest(fakeData,0,tParams{:});
            clustDistrib(p) = max(clustLength(fakeT));
        end
        critVal = prctile(clustDistrib,95); % 95 % of values are below this
   
        realLocs = bwlabel(realT);   %identify contiguous ones
        realLength = regionprops(realLocs, 'area');  %length of each span
        realLength = [ realLength.Area];
        corrIdx = find(realLength > critVal);
        corrT = ismember(realLocs,corrIdx);
        if deletePool
            delete(poolobj)
        else
        end

end

function rB = num2bits( aV, aN )
    tV = aV(:); % Make sure v is a row vector.
    tNV = length( tV ); % number of numbers to convert
    rB = zeros( tNV, aN );
    tP = aN - 1;
    rB( :, 1 ) = mod( tV, ( 2 ^ aN ) );
    for iP = 1:( aN - 1 )
        rB( :, iP+1 ) = mod( rB( :, iP ), ( 2 ^ tP ) );
        rB( :, iP ) = floor( rB( :, iP ) / ( 2 ^ tP ) );
        tP = tP - 1;
    end
    rB( :, end ) = floor( rB( :, end ) / ( 2 ^ tP ) );
    rB ( rB == 0 ) = -1;
end


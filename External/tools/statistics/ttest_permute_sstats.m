function [realH,realP,realT,corrH,critVal,supraTh,randDist]= ttest_permute_sstats( inData,maxPerms,testStat,timeMask,permMatrix,deletePool,tParams)
        % [realH,realP,corrH,critVal,clustDistrib]= ttest_permute_sstat( inData,maxPerms,timeMask,permMatrix,deletePool,tParams )
        %
        % 
        % input:
        %   inData: the difference waveform to test against zero, as a t x s matrix, where t is time and s is subjects.
        %   testStat: What kind of summary statistic to be used for clusters
        %             ['size']/'height'/'mass', (Cluster size/ Cluster height/ Cluster excess mass)
        %             for more info please refer to Pernet et al. 2015, Journal of Neuroscience Methods
        %             and Bullmore et al. 1999, IEEE transactions on medical imaging
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
        %   realH (1xt logical): significance of ordinary t-test
        %   realP (1xt double): p-values of ordinary t-test
        %   corrH (1xt logical): significance of corrected t-test
        %   critVal (scalar): critical number of significant time-points
        %                     required for surviving correction
        %   clustDistrib: distribution of significance run lengths for
        %                 permuted data.
         
% Latest modification Elham Barzegaran, 9/17/2018 -> to apply the cluster-based permutation test based on excess mass () of the cluster not its length:
% This allows short but strong differences to survive
%%  Set defaul values        
        % take out NaN data and count subjects
        notNaN = sum(isnan(inData),1)==0;
        inData = inData(:,notNaN);
        numSubs = size(inData,2);
        
        numPerms = 2 ^ numSubs;
        if ~exist('tParams','var') || isempty(tParams)
            tParams = {'dim',2,'alpha',0.05};
        else
        end
        if ~exist('deletePool','var') || isempty(deletePool)
            deletePool = false;
        else
        end
        if ~exist('permMatrix','var') || isempty(permMatrix)
            permMatrix = [];
        end
        if ~exist('timeMask','var') || isempty(timeMask)
            timeMask = 1:size(inData,1);
        else
        end
        if ~exist('maxPerms','var') || isempty(maxPerms)
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
        
        if ~exist('testStat','var') || isempty(testStat)
            testStat = {'size'};
        else
            if sum(strcmpi(testStat,{'size','height','mass'}))==0
                warning(['No summary statistic as ' testStat ' is implemented...']);
                warning('Continuing with cluster size as summary statistic.');
                testStat = {'size'};
            else
                testStat = lower(testStat);
            end
        end
        %%
        
        realData = inData(timeMask,:);
        realH = zeros(size(inData,1),1);
        realP = ones(size(inData,1),1);
        [realH(timeMask),realP(timeMask),~,realstats] = ttest(realData,0,tParams{:});
        realT = realstats.tstat;
        supraTh = abs(tinv(tParams{4}/2,realstats.df(1)));% suprathreshold for defining the clusters
        
        supraThm = abs(tinv(.2/2,realstats.df(1)));% suprathreshold for defining the clusters  
        supraThM = abs(tinv(min(realP)/2,realstats.df(1)));% suprathreshold for defining the clusters
        Thbins = supraThm:(supraThM-supraThm)/20:supraThM;        
        poolobj = gcp;
        
        parfor p = 1:maxPerms % all possible permutations
            fakeData = realData;
            fakeData = fakeData.*repmat(permMatrix(p,:),size(fakeData,1),1);
            [~,~,~,Stats] =  ttest(fakeData,0,tParams{:});
            %[clust_exessmass(p,:),clust_size(p,:),clust_height(p,:)] = arrayfun(@(x) clustStats(Stats.tstat,Thbins(x)),1:numel(Thbins));
            [clus1,clus2,clus3]= clustStats(Stats.tstat,supraTh);
            clust_mass(p) = max(clus1);clust_size(p) = max(clus2);clust_height(p) = max(clus3);
        end     
         
        realLocs = bwlabel(realH);   %identify contiguous ones
        realH2 = cat(find(max(size(realH))),0,realH,0);
        Clusters = [find(diff(realH2)==1) find(diff(realH2)==-1)-1];
        realCluster_size = arrayfun(@(x) Clusters(x,2)-Clusters(x,1)+1,1:size(Clusters,1)); % cluster size
        realCluster_mass = abs(arrayfun(@(x) sum(realstats.tstat(Clusters(x,1):Clusters(x,2))),1:size(Clusters,1)));% excess mass
        realCluster_mass = realCluster_mass - (realCluster_size*supraTh);
        realCluster_height = abs(arrayfun(@(x) max(realstats.tstat(Clusters(x,1):Clusters(x,2))),1:size(Clusters,1)));
        realCluster_height = realCluster_height - supraTh;
    
        % extract the significant clusters
        eval(['randDist = clust_' testStat ';']);
        eval(['realStat = realCluster_' testStat ';']);
        critVal = prctile(randDist,95); % 95 % of values are below this    
        corrIdx = find(realStat > critVal); % find clusters in real data bigger than 95% of null distribution
        corrH = ismember(realLocs,corrIdx);
 
        if false
            [Dist bins]= hist3([clust_size' clust_height'],'Edges',{1:max(max(clust_size),max(realCluster_size))+1,min(clust_height,2):.1:max(max(clust_height),max(realCluster_height))+1});
            figure,imagesc(Dist(end:-1:1,:));
            set(gca,'ytick',1:4:numel(bins{1}),'yticklabel',ceil(bins{1}(end:-4:1)'),'xtick',1:4:numel(bins{2}),'xticklabel',bins{2}(1:4:end)');
            xlabel('Cluster height');ylabel('cluster size')
        end
        if deletePool,delete(poolobj),end
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


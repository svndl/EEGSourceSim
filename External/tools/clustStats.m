function [Cluster_sstat,Cluster_size,Cluster_height] = clustStats(Stats,supraTh)
    % EB, 8/16/2018
    inData = abs(Stats)>supraTh;
    inData = cat(find(max(size(inData))),0,inData,0);
    Clusters = [find(diff(inData)==1) find(diff(inData)==-1)-1];
    
    Cluster_size = arrayfun(@(x) Clusters(x,2)-Clusters(x,1)+1,1:size(Clusters,1));
    
    Cluster_sstat = abs(arrayfun(@(x) sum(Stats(Clusters(x,1):Clusters(x,2))),1:size(Clusters,1)));% cluster mass according to (Bullmore et al., 1999)
    Cluster_sstat = Cluster_sstat - (Cluster_size*supraTh); % remove the cluster minimum height
    
    Cluster_height = abs(arrayfun(@(x) max(Stats(Clusters(x,1):Clusters(x,2))),1:size(Clusters,1)));% cluster height
    Cluster_height = Cluster_height - supraTh;
    
    if isempty(Cluster_sstat)
        Cluster_sstat   = 0;
        Cluster_size    = 0;
        Cluster_height  = 0;
    end
    
    % biggest cluster stats
    Cluster_size = max(Cluster_size);
    Cluster_sstat = max(Cluster_sstat);
    Cluster_height = max(Cluster_height);
end


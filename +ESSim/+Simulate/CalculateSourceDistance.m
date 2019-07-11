function [spat_dists, Euc_dist] = CalculateSourceDistance(MDATA,distanceType)
% Syntax: [spat_dists, Euc_dist] = CalculateSourceDistance(MDATA,distanceType)
% Description: This function calculates the distances between the sources using the
% distanceType method

%--------------------------------------------------------------------------
% INPUTS
    %   MDATA:  a structure which should have the following fields:
    %           - vertices: which is  a ns x 3 matrices indicating the source coordinates in mm
    %           - triangles: indicates the faces of the cortical meshe 
    %             (Necessary only if Geodesic distance is indicated in the input)
    %
    %   distanceType: indicates how to calculate the source distaces, [Euclidean]/ Geodesic 
    %   
%OUTPUTS
    %   spat_dists:  a ns x ns matrix indicating the sources in mm
%--------------------------------------------------------------------------
% Author: Elham Barzegaran
% Latest modification: 03.13.2018


%%
    Euc_dist = squareform(pdist(MDATA.vertices')) ;
    Euc_dist(round(length(Euc_dist)/2)+1:end,1:round(length(Euc_dist)/2))=Inf;
    Euc_dist(1:round(length(Euc_dist)/2),round(length(Euc_dist)/2)+1:end)=Inf;
    %---------------------EUCLIDEAN DISTANCE-------------------------------
    if strcmp(distanceType,'Euclidean') 
        spat_dists =  Euc_dist; 
        
    %---------------------GEODESIC DISTANCE--------------------------------   
    elseif strcmp(distanceType,'Geodesic')
        faces = MDATA.triangles'; vertex = MDATA.vertices';
        [c, ~] = tess_vertices_connectivity( struct( 'faces',faces + 1, 'vertices',vertex ) ); 
        % although using graph-based shortest path algorithm, we can estimate surface distances: But this will overestimate the real
        % distances...
        
        %---Using in-built matlab function:requires matlab 2015b ornewer---
        if exist('graph')==2
            G = c.*Euc_dist;G(isnan(G))=0;
            spat_dists = distances(graph(G));
            clear G c;
        else 
        %-----------Otherwise use Dijkstra function: very slow ------------
        % There are still some problem with this part....
            spat_dists = Euc_dist2;%ones(size(c))*max(Euc_dist(:));
            warning ('If the matlab version you are using is older than 2015b, Geodesic distances is not implmented. SUGGESTION: Use Euclidean distance instead');
            
            % this input gives an option to the user to switch to Euclidean if Geodesic is taking too much time

            inp = input('Do you want to continue with Euclidean distance?[Y]/N \n Enter Y to continue with Euclidean distance \n Enter N to continue with Geodesic distance\n','s');
            if strcmpi(inp,'N')
                error ('Geodesic distance is not implemented in this version of matlab')
            else
                spat_dists =  Euc_dist;
                return
            end
            %----------------------SLOW SHORTEST PATH----------------------
            hWait = waitbar(0,'Calculating Geodesic distances ... ');
            j = 0;% counter for waitbar
            for s = 1:length(c)
                if mod(s,20)==0, waitbar(j/length(c)); end
                if s<=MDATA.VertexLR(1) % Separate left and right hemispheres
                    sidx = 1:MDATA.VertexLR(1);
                else
                    sidx = 1+MDATA.VertexLR(1):sum(MDATA.VertexLR);
                end
                
                sidx = find(Euc_dist(s,:)<25);% or put a threshold
                if sum(sidx>s) % Apply diskstra algorithm to the hemisphere
                    spat_dists(s,sidx(sidx>=s)) = dijkstra(c(sidx,sidx),Euc_dist(sidx,sidx),find(sidx==s),find(sidx>=s),0);% complexity of this algorithm is O(|V^2|), where V is number of nodes
                end
                j = j+1;
 
            end
            spat_dists = min(spat_dists,spat_dists');
            close(hWait);
        end
    end 
end
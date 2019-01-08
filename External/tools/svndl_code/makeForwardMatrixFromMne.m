function [sol solFree] = makeForwardMatrixFromMne(fwd,srcSpace)
%function [sol] = makeForwardMatrixFromMne(fwd,srcSpace)
%
% fwd is the structure returned by mne_read_forward_solution
% srcSpace is the structure returned by mne_read_source_spaces
% sol is a matrix
%
%   $Log: makeForwardMatrixFromMne.m,v $
%   Revision 1.6  2010/07/22 17:01:53  SKI+ales
%   Fixed dumb typo
%
%   Revision 1.5  2010/07/21 23:57:37  SKI+ales
%   Fixed a nan bug
%
%   Revision 1.4  2009/11/12 22:14:48  ales
%   Fixed divide by zero error when too many sources omitted
%
%   Revision 1.3  2009/11/02 17:46:00  ales
%   Many bug fixes
%
%   Revision 1.2  2009/06/26 21:12:30  ales
%   mrSimScript does arbitrary waveforms
%   skeriDefaultSimParameters has updated help
%   makeDefaultCortexMorphMap returns an identity matrix for mapping an individual on to itself
%   makeForwardMatrixFromMne has some added help
%
%   Revision 1.1  2009/05/15 17:29:34  ales
%   Createsa a forward matrix from MNE FILES: a srcSpace and fwd
%

me='JMA:makeForwardMatrixFromMne';
FIFF=fiff_define_constants;

newDx = [];
oriLeft = [];
idx=1;
for iHi=fwd.src(1).vertno,
    newDx(idx) = find(srcSpace(1).vertno == iHi);
    oriLeft(idx,:) = fwd.src(1).nn(iHi,:);
    idx=idx+1;
    
end

newDxLeft =newDx;

newDx = [];
idx=1;
for iHi=fwd.src(2).vertno,
    newDx(idx) = find(srcSpace(2).vertno == iHi);
    oriRight(idx,:) = fwd.src(2).nn(iHi,:);
    idx=idx+1;
    
end

newDxRight =newDx;

totalDx = [newDxLeft (newDxRight+double(srcSpace(1).nuse))];




Gfree = fwd.sol.data;
nChan = size(Gfree,1);
Gfree = Gfree - repmat(mean(Gfree),nChan,1);

G = zeros(size(Gfree,1),size(Gfree,2)/3);
srcOri = [oriLeft; oriRight];

idx=1;
for i=1:3:length(Gfree);

    G(:,idx) = Gfree(:,i:i+2)*srcOri(idx,:)';

    idx=idx+1;
end

origSensPow = trace(G*G');

%This line is compeletely inscrutable
%Basically I am reducing the number of srcs by summing up a bunch of
%sources in a chunk (G*srcCor), then I am expanding that back to the
%orignal size of the matrix putting the summed source into it's place by
%remultiplying by srcCor


%[non G] = colnorm(G);

%newSensPow = trace(G*G');

%scaleFac = origSensPow/newSensPow;

%G = scaleFac*G; 


%trans = diag(sparse(inv.reginv))*inv.eigen_fields.data*inv.whitener*inv.proj;
%sol   = srcCov*inv.eigen_leads.data*trans;

% if inv.source_ori == FIFF.FIFFV_MNE_FREE_ORI
%     fprintf(1,'combining the current components...');
%     sol1 = zeros(size(sol,1)/3,size(sol,2));
%     for k = 1:size(sol,2)
%         sol1(:,k) = sqrt(mne_combine_xyz(sol(:,k)));
%     end
%     sol = sol1;
%    
% end

%sol = sol(3:3:end,:);

totalVerts = srcSpace(1).nuse + srcSpace(2).nuse;
sol1 = zeros(totalVerts,size(G,1));
solFree =zeros(3*totalVerts,size(G,1));

solFree([totalDx*3]-2,:) = Gfree(:,1:3:end)';
solFree([totalDx*3]-1,:) = Gfree(:,2:3:end)';
solFree(totalDx*3,:) = Gfree(:,3:3:end)';

sol1(totalDx,:) = G';


%The following averages values for the missing vertices


removedVerts = setdiff(1:srcSpace(1).nuse,newDxLeft);
for iGone = removedVerts,  
    [i j] = find(srcSpace(1).use_tris==srcSpace(1).vertno(iGone));
    neighbors = double(srcSpace(1).use_tris(i,:));
    neighbors = unique(neighbors(neighbors~=srcSpace(1).vertno(iGone)));
    neighbors = setdiff(neighbors,removedVerts);   
    [srcIdx] = find(ismember(srcSpace(1).vertno,neighbors));
    
    goodList = setdiff(srcIdx,removedVerts);
    
    if ~isempty(goodList)
        sol1(iGone,:) = mean(sol1(goodList,:),1);
        solFree(iGone*3-2,:) = mean(solFree(goodList*3-2,:),1);
        solFree(iGone*3-1,:) = mean(solFree(goodList*3-1,:),1);
        solFree(iGone*3,:)   = mean(solFree(goodList*3,:),1);
 
    end
   
    
    

end


removedVerts = setdiff(1:srcSpace(2).nuse,newDxRight);
for iGone = removedVerts,  
    [i j] = find(srcSpace(2).use_tris==srcSpace(2).vertno(iGone));
    neighbors = double(srcSpace(2).use_tris(i,:));
    neighbors = unique(neighbors(neighbors~=srcSpace(2).vertno(iGone)));
    neighbors = setdiff(neighbors,removedVerts);   
    [srcIdx] = find(ismember(srcSpace(2).vertno,neighbors));
    
    goodList = setdiff(srcIdx,removedVerts)+double(srcSpace(1).nuse);
    
    if ~isempty(goodList)
        sol1(iGone+double(srcSpace(1).nuse),:) = mean(sol1(goodList,:),1);
        solFree(iGone*3-2+3*double(srcSpace(1).nuse),:) = mean(solFree(goodList*3-2,:),1);
        solFree(iGone*3-1+3*double(srcSpace(1).nuse),:) = mean(solFree(goodList*3-1,:),1);
        solFree(iGone*3+3*double(srcSpace(1).nuse),:)   = mean(solFree(goodList*3,:),1);
    end

    


end

solFree(isnan(solFree(:))) = 0;
sol1(isnan(sol1(:))) = 0;

solFree = solFree';
sol = sol1';







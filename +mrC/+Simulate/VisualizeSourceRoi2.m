function [Fhandler,RoiList] = VisualizeSourceRoi2(subID,anatDir,RoiType,RoiIdx,direction,hemi,cmap,RoiLeg)
% gets the subject ID and anatomy folder and plots the
% ROIs on the subjects default cortex...
% Elham Barzegaran, 5.25.2018
%% default variables

if ~exist('anatDir','var')||isempty(anatDir)
    anatDir = getpref('mrCurrent','AnatomyFolder');
end

if ~exist('direction','var')||isempty(direction)
    direction = 'anterior';
end

if ~exist('hemi','var')
    hemi = 'B';
end

%%

% load default cortex
Path = fullfile(anatDir,subID,'Standard','meshes','defaultCortex.mat');
load(Path);
vertices = msh.data.vertices';
faces = (msh.data.triangles+1)';
% adjustements for visualization purpose
vertices = vertices(:,[1 3 2]);vertices(:,3)=200-vertices(:,3);

%% plot brain surface
if hemi=='L'
    Rind = round(length(vertices)/2)+1:length(vertices);
    [~,i1] = intersect(faces(:,1),Rind);[~,i2] = intersect(faces(:,2),Rind);[~,i3] = intersect(faces(:,3),Rind);
    I = ([i1; i2; i3]);
    faces(I,:)=[];
    [~,i1] = intersect(faces(:,1),Rind);[~,i2] = intersect(faces(:,2),Rind);[~,i3] = intersect(faces(:,3),Rind);
    I = ([i1; i2; i3]);
    faces(I,:)=[];
elseif hemi=='R'
    Lind = 1:round(length(vertices)/2);
    [~,i1] = intersect(faces(:,1),Lind);[~,i2] = intersect(faces(:,2),Lind);[~,i3] = intersect(faces(:,3),Lind);
    I = unique([i1; i2; i3]);
    faces(I,:)=[];
    [~,i1] = intersect(faces(:,1),Lind);[~,i2] = intersect(faces(:,2),Lind);[~,i3] = intersect(faces(:,3),Lind);
    I = unique([i1; i2; i3]);
    faces(I,:)=[];
end

%% plot ROIs on the brain
% RoiDir = fullfile(anatDir,subID,'Standard','meshes',[RoiType '_ROIs']); 
% [chunks,RoiList] = mrC.ChunkFromMesh(RoiDir,size(vertices,1),1);
% Fhandler= figure;


Rois = mrC.ROIs([],anatDir);
Rois = Rois.loadROIs(subID,anatDir);
Rois = Rois.getAtlasROIs(RoiType);
Rois = Rois.searchROIs('all',[],hemi);
chunks = Rois.ROI2mat(length(vertices));
RoiList = Rois.getFullNames('noatlashemi');

if ~exist('RoiIdx','var')||isempty(RoiIdx)
    RoiIdx = 1:size(chunks,2);
end
if ~exist('RoiLeg','var')% which ROI to present in legend: index of RoiIdx
    RoiLeg = 1:numel(RoiIdx);
end


if iscell(RoiIdx) && (ischar(RoiIdx{1}) || isstring(RoiIdx{1}))
    RoiIdx = cellfun(@(x) find(contains(RoiList,x)),RoiIdx);
end

if ~exist('cmap','var')
    cmap = distinguishable_colors(numel(RoiIdx),[.7 .7 .7]);
end
% cmap = cmap(randperm(numel(RoiIdx)),:);
isem = zeros(1,numel(RoiIdx));
pind = 1;
for i = 1:numel(RoiIdx)
    [RoiV, RoiF] = SurfSubsample(vertices, faces,find(chunks(:,RoiIdx(i))),'union');  
    C = cmap(i,:);
    if ~isempty(RoiF)
        hold on; 
        if ismember(i,RoiLeg)
            leg(pind) = patch('faces',RoiF,'vertices',RoiV,'edgecolor','none','facecolor','interp','facevertexcdata',repmat(C,size(RoiV,1),1),...
                 'Diffusestrength',.55,'AmbientStrength',.8,'specularstrength',.2,'FaceAlpha',.95,'facelighting','gouraud');
             pind = pind+1;
        else
            patch('faces',RoiF,'vertices',RoiV,'edgecolor','none','facecolor','interp','facevertexcdata',repmat(C,size(RoiV,1),1),...
                'Diffusestrength',.55,'AmbientStrength',.8,'specularstrength',.2,'FaceAlpha',.95,'facelighting','gouraud');
            %scatter3(RoiV(:,1),RoiV(:,2),RoiV(:,3),30,C,'filled');
        end
    else
        isem(i) = 1;
    end
end
%legend(leg,RoiList(RoiIdx(~isem))));



%%
[r,c] = find(chunks(:,RoiIdx));
ROIvers = sort(unique(r));
v1 = ismember(faces(:,1),ROIvers);
v2 = ismember(faces(:,2),ROIvers);
v3 = ismember(faces(:,3),ROIvers);
RemovVers = v1 .* v2 .*v3;

patch('faces',faces(~RemovVers,:),'vertices',vertices,'edgecolor','none','facecolor','interp','facevertexcdata',repmat([.72,.7,.7],size(vertices,1),1),...
     'Diffusestrength',.45,'AmbientStrength',.3,'specularstrength',.1,'FaceAlpha',.95,'facelighting','gouraud');

%colormap(cmap);

shading interp
lighting gouraud
lightangle(50,120)
lightangle(50,0)

switch direction
    case 'ventral'
        lightangle(-90,-90);
        view(-90,-90)
    case 'anterior'
        view (90,-10)
end

axis  off vis3d equal;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.24, .45, 0.65]);
if exist('leg','var') && ~isempty(leg)
    l = legend(leg,RoiList(RoiIdx(RoiLeg)),'location','west');
    set(gca,'fontsize',20)
    lp = get(l,'position');
    %set(l,'position',[0.0134    0.2423    0.1384    0.5504])
   % set(gca, 'position', get(gca,'position')+[0 0 lp(3:4)]);
end
end

function [nvertices, nfaces,vertIdx2] = SurfSubsample(vertices, faces,vertIdx,type)

% This function selects a subset of vertices and their corresponding faces indicated by vertIdx 
% INPUT: 
    % type: can be 'union' or 'intersect', the criteria for including faces
%-------------------------------------------------------------------------
if ~exist('type','var'),type = 'intersect';end
%-------------------------------------------------------------------------
vertIdx = sort(vertIdx);
I1 = find(ismember(faces(:,1),vertIdx));
I2 = find(ismember(faces(:,2),vertIdx));
I3 = find(ismember(faces(:,3),vertIdx));

if strcmp(type,'intersect')
    FI  = intersect(intersect(I1,I2),I3);
    nfaces = faces(FI,:);
    vertIdx2 = unique(nfaces(:));
    nvertices = vertices(vertIdx,:);
    fnew = zeros(size(nfaces));
    for i = 1:numel(vertIdx)
        fnew(nfaces==vertIdx(i)) = i;
    end
    nfaces=fnew;
    [~,vertIdx2] = intersect(vertIdx2,vertIdx);

elseif strcmp(type,'union')
    FI = union(union(I1,I2),I3);
    nfaces = faces(FI,:);
    vertIdx2 = unique(nfaces(:));
    nvertices = vertices(vertIdx2,:);
    fnew = zeros(size(nfaces));
    for i = 1:numel(vertIdx2)
        fnew(nfaces==vertIdx2(i)) = i;
    end
    nfaces=fnew;
    [~,vertIdx2] = intersect(vertIdx2,vertIdx);

else
    error('Criteria is not defined');
end

end
function Fhandler = VisualizeSourceData(subID,data,anatDir,cmap,direction,hemi)
% gets the subject ID and anatomy folder and the source data and plots the
% result on the subjects default cortex...
% Elham Barzegaran, 5.22.2018
%% default variables

if ~exist('anatDir','var') ||isempty(anatDir),
    anatDir = getpref('mrCurrent','AnatomyFolder');
end

if ~exist('direction','var') || isempty(direction),
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

if ~exist('cmap','var') || isempty(cmap)
    cmap = jmaColors('hotcortex');
end

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


Fhandler = patch('faces',faces,'vertices',vertices,'edgecolor','none','facecolor','interp','facevertexcdata',reshape((data),[numel(data) 1]),... 
     'Diffusestrength',.55,'AmbientStrength',.3,'specularstrength',.1,'facelighting','gouraud','FaceAlpha',.95);

colormap(cmap);

shading interp
%lighting flat
lightangle(50,120)
lightangle(50,0)

switch direction
    case 'ventral'
        lightangle(-90,-90);
        view(90,-90)
    case 'anterior'
        view (90,-10)
end

axis  off vis3d equal
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.24, .45, 0.65]);

end
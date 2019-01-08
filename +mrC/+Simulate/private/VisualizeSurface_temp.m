function VisualizeSurface(vertices,faces,FilePath,Filename,Hem,Data,ROI)

% Data should only contain possitive values for now
%%

if exist('Data','var')
    % remove Inf values
    Data(isinf(Data))=-1;
    Data((Data==-1))=max(Data(:));
    
    % normalize
    Data = ((Data-min(Data))/(max(Data)-min(Data)));
    Data = ((Data).^1)*63;
    %
    cmap = jmaColors('seedcortex');%jmaColors('Nebraska');
    %jet(round(max(Data(:)))+1);%hot(round(max(Data(:)))+1);
    Colors = cmap(round(Data)+1,:);
else
    Colors = repmat([.8,.8,.8],size(vertices,1),1);
end

switch Hem
    case 'L'
        Faces = faces(1:(size(faces,1))/2,:);
    case 'R'
        Faces = faces(((size(faces,1))/2)+1:end,:);
    case 'B'
        Faces = faces;
end
Brain =figure;

% Sources
if exist('ROI','var')
    scatter3(vertices(ROI,1),vertices(ROI,2),vertices(ROI,3),100,'c','filled');
end

% Brain surface
patch('faces',Faces,'vertices',vertices,'edgecolor','none','facecolor','interp','facevertexcdata',Colors,...
     'Diffusestrength',.45,'AmbientStrength',.3,'specularstrength',.1,'FaceAlpha',.95);
 
shading interp
lighting flat
lightangle(50,120)
%lightangle(50,240)
lightangle(50,0)
view(90,0)
axis  off vis3d equal
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.24, .45, 0.65]);
view(90,0);
im{1} = getframe(Brain);
for i = 1:179
    view(90-(i*2),0);
    im{i+1} = getframe(Brain);
end

%% create the video writer
 writerObj = VideoWriter(fullfile(FilePath,[Filename '.avi']));
 writerObj.FrameRate = 20;
 writerObj.Quality = 100;
 % open the video writer
 open(writerObj);

 % write the frames to the video
 for u=2:length(im)
     %frame = im2frame(im{u});
     writeVideo(writerObj, im{u});
 end

 % close the writer object
 close(writerObj);
end
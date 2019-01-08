function [ctx] = readDefaultCortex(subjId)
%[ctx] = readDefaultCortex(subjId)


anatDir = getpref('mrCurrent','AnatomyFolder');

ctxFilename=fullfile(anatDir,subjId,'Standard','meshes','defaultCortex');

load(ctxFilename);

ctx.faces = (msh.data.triangles+1)';

ctx.vertices(:,1) = msh.data.vertices(3,:)-128;
ctx.vertices(:,2) = -msh.data.vertices(1,:)+128;
ctx.vertices(:,3) = -msh.data.vertices(2,:)+128;
ctx.vertices = ctx.vertices/1000;

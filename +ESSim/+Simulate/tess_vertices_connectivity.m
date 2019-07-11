function [C, VertConn] = tess_vertices_connectivity(FV,VERBOSE)
% Computes vertices connectivity fast.
%
% USAGE:  [C, VertConn] = tess_vertices_connectivity(FV,VERBOSE)
% 
% INPUT:
%    - FV      : The standard matlab structure for faces and vertices,
%                where FV.faces is m x 3, and FV.vertices is n x 3.
%    - VERBOSE : (optional) used for compatibility with an earlier function
% OUTPUT:
%    - VertConn: vector of cells, one-to-one with each row of FV.vertices.
%                VertConn{i} returns a row vector of the vertex numbers (rows in FV.vertices) that
%                are connected by faces to the ith vertex row in FV.vertices.
%                Thus if we want to 'swell' the region around a vertex, VertConn{i} gives us the 
%                vertex numbers of those vertices that are adjacent.
%    - C       : Connectivity sparse matrix with dimension nVertices x nVertices. It has 1
%                at (i,j) when vertices i and j are connected

% @=============================================================================
% This software is part of The BrainStorm Toolbox
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2008 BrainStorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm licence" at command prompt.
% =============================================================================@
%
% Authors: Anand Joshi, Dimitrios Pantazis, November 2007 
% ----------------------------- Script History ---------------------------------
% DP  01-Nove-2007  Creation
% ------------------------------------------------------------------------------


%make connectivity matric C
rowno = double([FV.faces(:,1);FV.faces(:,1);FV.faces(:,2);FV.faces(:,2);FV.faces(:,3);FV.faces(:,3)]);
colno = double([FV.faces(:,2);FV.faces(:,3);FV.faces(:,1);FV.faces(:,3);FV.faces(:,1);FV.faces(:,2)]);
data=ones(size(rowno));
C = logical(sparse(rowno,colno,data));

if (nargout >= 2)
    %find nonzero elements and initialize
    [rows,cols] = find(C);
    d = find(diff([cols' cols(end)+1]));
    nVertices = size(FV.vertices,1);
    VertConn = cell(nVertices,1); % one empty cell per vertex

    % If there are vertices without any connections, do a slow loop
    if length(d)~=nVertices
        for i = 1:length(cols)
            VertConn{cols(i)} = [VertConn{cols(i)} rows(i)];
        end
    % Else: fast calculation of vertices connectivity (if no isolated vertices)
    else
        VertConn{1} = rows(1:d(1))';
        for i = 2:nVertices
            VertConn{i} = rows((d(i-1)+1):d(i))';
        end
    end
else
    VertConn = [];
end



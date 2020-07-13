function [pmesh_new, tmesh_new, bedges_new] = TriangularMeshRefinement2D(pmesh, tmesh, bedges)

%% Function and parameter description
% Helper function for unstructured grid refinement
% This function takes a triangular mesh and refines it
%
% The algorithm divides each triangle into 4 new triangles
% Every boundary edge is splitted apart into 2 new edges
%
%  1-------3           1---6---3
%   \     /             \ / \ /   
%    \   /     ---->     4---5
%     \ /                 \ /
%      2                   2
%      
%   ______________             _______________
%  |\            /|           |\     / \     /|
%  |  \        /  |           |  \ /_____\ /  | 
%  |    \    /    |           |  /|\     /|\  |
%  |      \/      |  ----- >  |/  |  \_/  |  \|
%  |      /\      |           |\  |  / \  |  /|
%  |    /    \    |           |  \|/_____\|/  |
%  |  /        \  |           |  / \     / \  |
%  |/____________\|           |/_____\_/_____\|


% Returns:
% pmesh_new -> new point mesh with unique coordinates
% tmesh_new -> new triangle mesh with unique coordinates

% Input args:
% pmesh, tmesh -> input mesh that shall be divided

%% Function logic

% Convert to double for using triangulation
tmesh = double(tmesh); 

% Validate input arguments

if (nargin ~= 3)
    error('Wrong number of arguments, put in pmesh, tmesh and bedges');
end

if (size(pmesh, 2) < 2) 
    error('Vertices should contain at least 2 columns');
elseif (size(tmesh, 2) < 3)
    error('Triangles should contain at least 3 columns');
elseif (size(bedges, 2) < 3)
    error('BEdges should contain 3 columns, third column is a value that gets copied on new edge');
end



numP = size(pmesh,1); % number of points
numT = size(tmesh,1); % number of triangles

Edges = zeros(3 * numT, 2);  % every triangle has 3 edges
edgeMid = false(3 * numT,1);
e = [0 0;0 0;0 0];
r = [0 0 0];

% generate triMesh struct with triangulation
triMesh = triangulation(tmesh,pmesh);
Edges = edges(triMesh);

midPoints = (pmesh(Edges(:,1),:) + pmesh(Edges(:,2),:))/2; % mid points
numMidPoints = size(midPoints, 1);

numNewP = numP + numMidPoints; % number of new points
numNewT = 4 * numT;            % number of new triangles

tmesh_new = zeros(numNewT, 3); 
pmesh_new = zeros(numNewP, 2);

pmesh_new(1:numP,:) = pmesh;
pmesh_new((numP+1):(numP+numMidPoints),:) = midPoints;

rowHelper = [1 1 1];
oldPoints = [1 1 1];

% Split every triangle into 4 new triangles
for n = 1 : numT
    aface = tmesh(n,:);
    
    %find and sort edge in a triangle
    if aface(1) < aface(2)
        e(1,1) = aface(1);
        e(1,2) = aface(2);
    else
        e(1,2) = aface(1);
        e(1,1) = aface(2);
    end
    
    if aface(2) < aface(3)
        e(2,1) = aface(2);
        e(2,2) = aface(3);
    else
        e(2,2) = aface(2);
        e(2,1) = aface(3);
    end
    
    if aface(1) < aface(3)
        e(3,1) = aface(1);
        e(3,2) = aface(3);
    else
        e(3,2) = aface(1);
        e(3,1) = aface(3);
    end
    
    rowHelper(1) = findThisRow(Edges, e(1,:));
    rowHelper(2) = findThisRow(Edges, e(2,:));
    rowHelper(3) = findThisRow(Edges, e(3,:));
    
    oldPoints(1:3) = tmesh(n,:);
    
    tmesh_new(4*(n-1)+1,1) = oldPoints(1);
    tmesh_new(4*(n-1)+1,2) = numP + rowHelper(1);
    tmesh_new(4*(n-1)+1,3) = numP + rowHelper(3);
    
    tmesh_new(4*(n-1)+2,1) = numP + rowHelper(1);
    tmesh_new(4*(n-1)+2,2) = oldPoints(2);
    tmesh_new(4*(n-1)+2,3) = numP + rowHelper(2);
    
    tmesh_new(4*(n-1)+3,1) = numP + rowHelper(2);
    tmesh_new(4*(n-1)+3,2) = oldPoints(3);
    tmesh_new(4*(n-1)+3,3) = numP + rowHelper(3);
    
    tmesh_new(4*(n-1)+4,1) = numP + rowHelper(1);
    tmesh_new(4*(n-1)+4,2) = numP + rowHelper(2);
    tmesh_new(4*(n-1)+4,3) = numP + rowHelper(3);
end


%% Adding new boundary edges

% Create new Midpoint and find that point in newly generated pmesh
% Then add the old bcondition value to the new and old edge

%% TODO: HERE IS ERROR !!!!!!

numEdges = size(bedges,1) ;
bedges_new = zeros(numEdges, 3);

          
for i=1:numEdges      
    p1 = bedges(i,1);
    p2 = bedges(i,2);
    b_val = bedges(i,3); % the value of boundary
    
    x1 = pmesh(p1,1);
    x2 = pmesh(p2,1);
    
    y1 = pmesh(p1,2);
    y2 = pmesh(p2,2);
    
    mx = (x1+x2)/2;
    my = (y1+y2)/2;
    
    midPoint = [mx my];
    
    index = find(ismember(pmesh_new, midPoint, 'rows'));
    
    bedges_new(i,:) =        [p1, index, b_val];
    bedges_new(i+numEdges,:) = [index, p2, b_val];
    
end % for


end

%% Helper function to find a specific row in a matrix

function idx = findThisRow(matrixToSearch, rowToFind)

% find index of specific row in a matrix

numRows = size(matrixToSearch, 1);
idx = 1;

indexFound = false;
count = 1;

while (count <= numRows) && (indexFound == false)
    if (matrixToSearch(count,1) == rowToFind(1)) && (matrixToSearch(count,2) == rowToFind(2))
        idx = count;
        indexFound = true;
    end
    count = count + 1;
end

end





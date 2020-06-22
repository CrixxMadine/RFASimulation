function [pmesh_new, tmesh_new, bedges_new] = TriangularMeshRefinement2D(pmesh, tmesh, test, bedges)

%% Function and parameter description
% Helper function for unstructured grid refinement
% This function takes a triangular mesh and refines it
% The algorithm divides each triangle into 4 new triangles
%  ________             ________
%  \      /             \  /\  /   
%   \    /     ---->     \/__\/
%    \  /                 \  /
%     \/                   \/

%   ________________             _______________
%  |\              /|           |\     / \     /|
%  |  \          /  |           |  \ /_____\ /  | 
%  |     \    /     |           |  /|\     /|\  |
%  |       \/       |  ----- >  |/  |  \_/  |  \|
%  |       /\       |           |\  |  / \  |  /|
%  |     /    \     |           |  \|/_____\|/  |
%  |  /          \  |           |  / \     / \  |
%  |/______________\|           |/_____\_/_____\|


% Returns:
% pmesh_new -> new point mesh with unique coordinates
% tmesh_new -> new triangle mesh with unique coordinates

% Input args:
% pmesh, tmesh -> input mesh that shall be divided

%% Function logic



% function [W, G] = meshSubdivision(V, F)
if (test == 0)
extra = zeros(size(pmesh,1),1);
V = [pmesh, extra];
else 
V = pmesh;
end

F = double(tmesh);

%Output faces G are unique.
% V is n by 3 array of input vertices  
% F is m by 3 array of input faces
% W is array of output vertices 
% G is array of output faces
%
%      v1                   v1
%     / \                  / \
%    /   \      ->        a---c
%   /     \              / \ / \
% v2 ----- v3          v2-- b --v3
%
% start with checking input arguments

% if nargin ~= 2
%     error('wrong number of arguments');
% end

if (size(V, 2) ~= 3) 
    error('vertices should contain 3 columns');
elseif (size(F, 2) ~= 3)
    error('faces should contain 3 columns');
end
nV = size(V,1); %number of vertices
nF = size(F,1); %number of faces
Edges = zeros(3 * nF, 2); %assuming 3 edges on each face
edgeMid = false(3 * nF,1);
e = [0 0;0 0;0 0];
r = [0 0 0];
triMesh = triangulation(F,V);
Edges = edges(triMesh);
MidPoints = (V(Edges(:,1),:) + V(Edges(:,2),:))/2; %mid points
nMidPoints = size(MidPoints, 1);
nW = nV + nMidPoints; %number of new vertices
nG = 4 * nF; %number of new faces
tmesh_new = zeros(nG, 3); %new faces array
pmesh_new = zeros(nW, 3); %new vertices array
pmesh_new(1:nV,:) = V;
pmesh_new((nV+1):(nV+nMidPoints),:) = MidPoints;
k = [1 1 1];
ind = [1 1 1];
for n = 1 : nF
    aface = F(n,:);
    %find and sort edged in a face
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
    
    k(1) = findaRow(Edges, e(1,:));
    k(2) = findaRow(Edges, e(2,:));
    k(3) = findaRow(Edges, e(3,:));
    
    ind(1:3) = F(n,:);
    
    tmesh_new(4*(n-1)+1,1) = ind(1);
    tmesh_new(4*(n-1)+1,2) = nV + k(1);
    tmesh_new(4*(n-1)+1,3) = nV + k(3);
    
    tmesh_new(4*(n-1)+2,1) = nV + k(1);
    tmesh_new(4*(n-1)+2,2) = ind(2);
    tmesh_new(4*(n-1)+2,3) = nV + k(2);
    
    tmesh_new(4*(n-1)+3,1) = nV + k(2);
    tmesh_new(4*(n-1)+3,2) = ind(3);
    tmesh_new(4*(n-1)+3,3) = nV + k(3);
    
    tmesh_new(4*(n-1)+4,1) = nV + k(1);
    tmesh_new(4*(n-1)+4,2) = nV + k(2);
    tmesh_new(4*(n-1)+4,3) = nV + k(3);
end


%% TRY Adding new bedges -

% Create new Midpoint and find that point in new pmesh
% Then add old value to both edges

if (test == 1)

nEdges = size(bedges,1) ;
 
bedges_new = zeros(nEdges, 3);

          
for i=1:nEdges      
    p1 = bedges(i,1);
    p2 = bedges(i,2);
    b_val = bedges(i,3); % the value of boundary
    
    x1 = pmesh(p1,1);
    x2 = pmesh(p2,1);
    
    y1 = pmesh(p1,2);
    y2 = pmesh(p2,2);
    
    mx = (x1+x2)/2;
    my = (y1+y2)/2;
    
    midPoint = [mx my 0];
    
    index = find(ismember(pmesh_new, midPoint, 'rows'));
    
    bedges_new(i,:) =        [p1, index, b_val];
    bedges_new(i+nEdges,:) = [index, p2, b_val];
    
end % for
end % if test == 1

end

%% Helper

function k = findaRow(arr1, arr2)
%find arr2 that is row vector in another vector arr1
rows = size(arr1, 1);
k = 1;
indexFound = false;
n = 1;
while (n <= rows) && (indexFound == false)
    if (arr1(n,1) == arr2(1)) && (arr1(n,2) == arr2(2))
        k = n;
        indexFound = true;
    end
    n = n + 1;
end

end





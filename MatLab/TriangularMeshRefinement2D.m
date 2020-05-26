function [refinedMesh] = TriangularMeshRefinement2D(inputMesh)

% Helper function for structured grid refinement
% This function takes a triangular mesh and refines it uniformly
% The algorithm divides each triangle uniformly into 4 triangles
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


% TODO

refinedMesh = inputMesh;


end


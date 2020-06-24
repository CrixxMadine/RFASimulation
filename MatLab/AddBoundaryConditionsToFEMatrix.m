function [Ah_bound, fh_bound] = AddBoundaryConditionsToFEMatrix(Ah, fh, pmesh, tmesh, bmesh)

%% Function summary and arguments description

% Adds the specific boundary conditions for the current system of equations

% returns
% Ah_bound := FE-matrix with new boundary conditions
% fh_bound := right hand side with new boundary conditions

% Input args:
% Ah := FE-matrix (sum of stiffnessmatrix and massmatrix) 
% fh := right hand side of the system of equations
% pmesh  := point matrix of triangulation
% tmesh  := triangle matrix of triangulation

% bmesh  := boundary edge matrix with boundary points
%   -> First and second rows : pair of points of the boundary edge
%   -> Third row  :  boundary type (1 = dirichlet, 2 = neumann)
%   -> Fourth row :  right hand side value of the boundary condition


boundaryNodes = unique([bmesh(1,:), bmesh(2,:)]);
boundaryNodesNumber = size(boundaryNodes,2);

totalNodesNumber = size(pmesh,2);
innerNodesNumber = totalNodesNumber - boundaryNodesNumber;

innerNodes = setdiff(1:totalNodesNumber, boundaryNodes);

% Inlcude boundary condition

% Sorting bmesh
[~,idx] = sort(bmesh(1,:));   % sort just the first row
sorted_bmesh = bmesh(:,idx);  % sort the whole matrix using the sort indices

% First: Only dirichlet
dirichletValues = zeros(totalNodesNumber,1);
%isRealDirich = zeros(totalNodesNumber,1);

for i=1:boundaryNodesNumber 
    type = sorted_bmesh(3,i);
    
    if (type == 1) % Dirichlet
        node = boundaryNodes(i);
        dirichletValues(node) = sorted_bmesh(4,i); 
        
        % Testing with flag
        %isRealDirich(node) = 1;
    
        % TODO TRYING
%         Ah(node,node) = 1;
%         fh(node) = sorted_bmesh(4,i);
%         
%         Ah(node, node+1) = 0; 
%         Ah(node, node+2) = 0; 
%         
%         Ah(node+1, node) = 0; 
%         Ah(node+2, node) = 0; 
%         
%         fh(node+1) = fh(node+1) - Ah(node+1,node) * sorted_bmesh(4,i);
%         fh(node+2) = fh(node+2) - Ah(node+2,node) * sorted_bmesh(4,i);
        % END TODO
        
    elseif (type == 2) % Neumann
        % TODO This is wrong! -> Does Dirichlet conditions with value = 0
%         node = boundaryNodes(i);
%         dirichletValues(node) = 37 + 273.15;
%         isRealDirich(node) = 1;
        
        test = 0;
        
    else 
        % TODO
        
    end
end

% Ah_bound = Ah;
% fh_bound = fh;

% Transform matrix and load vector with dirichlet conditions
fh_bound = fh - Ah * dirichletValues;
fh_bound(boundaryNodes) = dirichletValues(boundaryNodes);

Ah_bound = Ah;

Ah_bound(boundaryNodes, 1:totalNodesNumber) = 0;
Ah_bound(1:totalNodesNumber, boundaryNodes) = 0;
for s=1:boundaryNodesNumber
   Ah_bound(boundaryNodes(s),boundaryNodes(s)) = 1;
end


%% TODO TESTING!!
% fh_bound = fh - Ah * dirichletValues;


% Test with for
% for j=1:boundaryNodesNumber
%     node = boundaryNodes(j);
%     
%     if (isRealDirich(node) == 1)
%         fh_bound(boundaryNodes) = dirichletValues(boundaryNodes);
%     end
% end
%
% Ah_bound = Ah;

% for zz=1:boundaryNodesNumber
%     node = boundaryNodes(zz);
%     
%     if (isRealDirich(node) == 1)
%         Ah_bound(node, 1:3) = 0;
%         Ah_bound(1:3, node) = 0;
%         Ah_bound(node,node) = 1;
%     end
% end

%Ah_bound(boundaryNodes, 1:totalNodesNumber) = 0;
%Ah_bound(1:totalNodesNumber, boundaryNodes) = 0;
% 
% for s=1:boundaryNodesNumber
%    %bNode = boundaryNodes(s);
%    
%    %if  bmesh(bNode)
%        Ah_bound(boundaryNodes(s),boundaryNodes(s)) = 1;
%    %end
% end



end


function [Ah_bound, fh_bound] = AddBoundaryConditionsToFEMatrix(Ah, fh, pmesh, bmesh)

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


boundaryNodes = unique([bmesh(:,1), bmesh(:,2)]);
boundaryNodesNumber = size(boundaryNodes,1);

totalNodesNumber = size(pmesh,1);

% Not important by now -> TODO
%innerNodesNumber = totalNodesNumber - boundaryNodesNumber;
%innerNodes = setdiff(1:totalNodesNumber, boundaryNodes);

% Inlcude boundary condition

% Sorting bmesh
[~,idx] = sort(bmesh(:,1));   % sort just the first column
sorted_bmesh = bmesh(idx,:);  % sort the whole matrix using the sort indices


%[~,idx] = sort(bmesh(:,2));
%sorted_bmesh_2 = bmesh(idx,:);


allNodes = unique([bmesh(:,1), bmesh(:,2)]);
left_Nodes = unique(bmesh(:,1));
right_Nodes = unique(bmesh(:,2));

notInLeft = setdiff(allNodes, left_Nodes);
notInRight = setdiff(allNodes, right_Nodes);


hereIsDirich = find(bmesh(:,3) == 1);

theseAreDirichNodes = unique([bmesh(hereIsDirich,1), bmesh(hereIsDirich,2)]);

% First: Only dirichlet
dirichletValues = zeros(totalNodesNumber,1);
isRealDirich = zeros(totalNodesNumber,1);
neumannValues = zeros(totalNodesNumber,1);

for i=1:boundaryNodesNumber 
    
    % Dirich Values
    node = boundaryNodes(i);
    if (ismember(node, theseAreDirichNodes))
        
         % THIS IS PROBLEM -> sets wrong number to dirich
         
         isRealDirich(node) = 1;
         

         
         % Is correct dirich value for node
         if (sorted_bmesh(i,3) == 1)
            dirichletValues(node) = sorted_bmesh(i,4);  
         
         else
             row = find(sorted_bmesh(:,1) == node);
             if ((sorted_bmesh(row,3) == 1))
              dirichletValues(node) = sorted_bmesh(row,4);
             else
         
         % Need to find correct dirich vale
            row = find(sorted_bmesh(:,2) == node);
            dirichletValues(node) = sorted_bmesh(row,4);
             end
         end
         
    
    % No dirich Value
    else
        
        isRealDirich(node) = 2;
        type = sorted_bmesh(i,3);
        
        if (type == 2) % Neumann           
            
            neumannValues(node) = sorted_bmesh(i,4);  
            
        else        
            % TODO        
        end
    %end
%     type = sorted_bmesh(i,3);
%     node = boundaryNodes(i);
%     
%     if (type == 1) % Dirichlet
% 
%         dirichletValues(node) = sorted_bmesh(i,4);       
%         isRealDirich(node) = 1;
   
        

        
    end
end


%% New Implementation -> Combine dirich and neumann values

% Calculate right hand side with substraction for dirichlet values
% (dirichvalues is zero, if it is not a dirichlet value)

Ah_bound = Ah;
fh_bound = fh - Ah * dirichletValues;

% Add all boundary conditions
for j=1:boundaryNodesNumber
    node = boundaryNodes(j);
    
    if (isRealDirich(node) == 1) % Dirichlet value
        Ah_bound(node, 1:totalNodesNumber) = 0;
        Ah_bound(1:totalNodesNumber, node) = 0;
        
        Ah_bound(node,node) = 1;
        
        fh_bound(node) = dirichletValues(node);
        
    elseif (isRealDirich(node) == 2) % Neumann value
         fh_bound(node) = fh(node) + neumannValues(node);   
         
    else 
        % TODO
    end
end

%% Old implementation  -> Assuming only dirichlet conditions
% Transform matrix and load vector with dirichlet conditions
%fh_bound = fh - Ah * dirichletValues;
% fh_bound(boundaryNodes) = dirichletValues(boundaryNodes);
% 
% Ah_bound = Ah;
% 
% Ah_bound(boundaryNodes, 1:totalNodesNumber) = 0;
% Ah_bound(1:totalNodesNumber, boundaryNodes) = 0;
% for s=1:boundaryNodesNumber
%    Ah_bound(boundaryNodes(s),boundaryNodes(s)) = 1;
% end



end


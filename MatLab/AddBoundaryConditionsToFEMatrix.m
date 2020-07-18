function [Ah_bound, fh_bound] = AddBoundaryConditionsToFEMatrix(Ah, fh, pmesh, bmesh)

%% Function summary and arguments description

% Adds the specific boundary conditions for the current system of equations
%
% returns
%  Ah_bound := FE-matrix with new boundary conditions
%  fh_bound := right hand side with new boundary conditions
%
% Input args:
%  Ah := FE-matrix (sum of stiffnessmatrix and massmatrix) 
%  fh := right hand side of the system of equations
%  pmesh  := point matrix of triangulation
%  tmesh  := triangle matrix of triangulation
%
%  bmesh  := boundary edge matrix with boundary points
%   -> First and second rows : pair of points of the boundary edge
%   -> Third row  :  boundary type (1 = dirichlet, 2 = neumann)
%   -> Fourth row :  right hand side value of the boundary condition


%% Implementation


% 1.) -> Find all boundary information in bmesh

boundaryNodes = unique([bmesh(:,1), bmesh(:,2)]);
boundaryNodesNumber = size(boundaryNodes,1);

totalNodesNumber = size(pmesh,1);


% Sorting bmesh
[~,idx] = sort(bmesh(:,1));   % sort just the first column
sorted_bmesh = bmesh(idx,:);  % sort the whole matrix using the sort indices


hereIsDirich = find(bmesh(:,3) == 1);
theseAreDirichNodes = unique([bmesh(hereIsDirich,1), bmesh(hereIsDirich,2)]);


dirichletValues = zeros(totalNodesNumber,1);
isRealDirich = zeros(totalNodesNumber,1);
neumannValues = zeros(totalNodesNumber,1);

% 2.) -> Loop through boundary edges and extract boundary values

for i=1:boundaryNodesNumber 
    
    node = boundaryNodes(i);
    
    % Case node is a dirichlet node
    if (ismember(node, theseAreDirichNodes))
               
         isRealDirich(node) = 1;
         
         row1 = find(sorted_bmesh(:,1) == node);
         row2 = find(sorted_bmesh(:,2) == node);
         row = [row1 row2];
             
         % Need to find correct dirich value    
         if ((sorted_bmesh(row(1),3) == 1))
             
            dirichletValues(node) = sorted_bmesh(row(1),4);

         elseif ((sorted_bmesh(row(2),3) == 1))
                      
            dirichletValues(node) = sorted_bmesh(row(2),4);

         end
         
    
    % Case node is not a dirichlet value
    else % if(ismember(node, dirich))
        
        isRealDirich(node) = 2;
        type = sorted_bmesh(i,3);

        if (type == 2) % Neumann           

            neumannValues(node) = sorted_bmesh(i,4);  

        else
            % Reserved for Robin conditions
        end
        
    end
end


% 3.) -> Add boundary conditions to the FEM system of equations

% Calculate right hand side with substraction for dirichlet values
% (dirichValues is set to zero, if it is not a dirichlet value)

Ah_bound = Ah;
fh_bound = fh - Ah * dirichletValues;

% Add all boundary conditions
for j=1:boundaryNodesNumber
    node = boundaryNodes(j);
    
    if (isRealDirich(node) == 1) % Dirichlet value
        
        Ah_bound(1:totalNodesNumber, node) = 0;

        Ah_bound(node, 1:totalNodesNumber) = 0;

        Ah_bound(node,node) = 1;
        
        fh_bound(node) = dirichletValues(node);
        
    elseif (isRealDirich(node) == 2) % Neumann value
        
         fh_bound(node) = fh(node) + neumannValues(node);   
         
    else 
        % Reserved for Robin conditions
    end
end



end


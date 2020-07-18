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


%% Implementation

boundaryNodes = unique([bmesh(:,1), bmesh(:,2)]);
boundaryNodesNumber = size(boundaryNodes,1);

totalNodesNumber = size(pmesh,1);


% Sorting bmesh
[~,idx] = sort(bmesh(:,1));   % sort just the first column
sorted_bmesh = bmesh(idx,:);  % sort the whole matrix using the sort indices


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
        % if (sorted_bmesh(i,3) == 1)
         %   dirichletValues(node) = sorted_bmesh(i,4);  
         
        % else
        
            % TODO -> this not right ...
             row1 = find(sorted_bmesh(:,1) == node);
             row2 = find(sorted_bmesh(:,2) == node);
             row = [row1 row2];
             
             testVal = sorted_bmesh(row,3);
             if (testVal(1) ~= testVal(2))
                stopHere = 0; 
                
                node1=sorted_bmesh(row(1),1);
                node2=sorted_bmesh(row(2),1);
                node3=sorted_bmesh(row(1),2);
                node4=sorted_bmesh(row(2),2);
                
                testestest = unique([node1 node2 node3 node4]);
                %isRealDirich(testestest) = 1;
             end
             
             if ((sorted_bmesh(row(1),3) == 1))
              dirichletValues(node) = sorted_bmesh(row(1),4);
              
%               if (testVal(1) ~= testVal(2))
%               % Testing
%               dirichletValues(testestest) = sorted_bmesh(row(1),4);
%               end
              
             elseif ((sorted_bmesh(row(2),3) == 1))
                          % Need to find correct dirich vale
            %row = find(sorted_bmesh(:,2) == node);
          dirichletValues(node) = sorted_bmesh(row(2),4);
          
          
%           if (testVal(1) ~= testVal(2))
%                         % Testing
%           dirichletValues(testestest) = sorted_bmesh(row(2),4);
%          end
          
             else 
                 error('Whaaaaa');

           
             end
       %  end
         
    
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
        
        Ah_bound(1:totalNodesNumber, node) = 0;

        Ah_bound(node, 1:totalNodesNumber) = 0;

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


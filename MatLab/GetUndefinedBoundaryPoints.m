function undefinedPoints = GetUndefinedBoundaryPoints(bmesh)

%% Function and parameter description

% Find problematic boundary nodes in grid
% 
%  Problematic node happens if there are different boundary conditions ...
%  colliding at one point e.g. the electrodes on the probe
%
%  Since the dirichlet value has higher priority and therefore set the its ... 
%  value, the next neighbour on the neumann or robin side is calculated wrong
% 
%  This functions purpose is to detect these problematic nodes
%
% Input arg:
%  bmesh, where each row represents a boundary edge
%  First and second column are the bnodes, third column is boundary type
%  Boundary types: 1 -> Dirichlet, 2 -> Neumann, 3 -> Robin
%
% Returns:
%  undefinedPoints: -> n x 3 matrix with the problematic node information
%  First column contains the problematic point in pmesh
%  Second and third column are the neigbour nodes of the undefined node


%% Implementation 

% Step 1.) ->  Find dirichlet nodes with non dirichlet neighbour if any

rowsDirichCondition = find(bmesh(:,3) == 1);
allDirichNodes = unique([bmesh(rowsDirichCondition,1), bmesh(rowsDirichCondition,2)]);

rowsNeumannCondition = find(bmesh(:,3) ~= 1);
allNeumannNodes = unique([bmesh(rowsNeumannCondition,1), bmesh(rowsNeumannCondition,2)]);

% Dirich node with neighbour that has other boundary condition
doubleDefinedNodes = intersect(allDirichNodes, allNeumannNodes);

numDoubleDefinedPoints = size(doubleDefinedNodes,1);

% Case no problem dirich node (e.g. only one type of boundary conditions)
if (numDoubleDefinedPoints == 0) 
    undefinedPoints = zeros(0,0);
    return;
end 


% Step 2.) -> Find all problematic neumann oder robin neighbours

problematicNodes(numDoubleDefinedPoints) = 0;

for count=1:numDoubleDefinedPoints
    
    corner = doubleDefinedNodes(count);
    
    rowHelper = [find(bmesh(:,1) == corner) find(bmesh(:,2) == corner)];    
    rowOfDirichCorner = unique(rowHelper);

    pointsNextToDirichCorner = [bmesh(rowOfDirichCorner,1) bmesh(rowOfDirichCorner,2)];

    cornerNeighbours = setdiff(pointsNextToDirichCorner, allDirichNodes);

    % On very coarse grid, this node might be defined     
    if (isempty(cornerNeighbours))
        cornerNeighbours = 0;
    end
    
    problematicNodes(count) = cornerNeighbours;
    
end

% On coarse grid it is possible, there is only one problem node ...
% in between two dirichlit nodes, e.g. between the electrodes
% Given that case, number of problematic nodes reduces
temp = unique(problematicNodes);
uniqueProblemNodes = temp(temp~=0);


% 3.) Find neighbours of problematic nodes and set values to return object

numProblemPoints = length(uniqueProblemNodes);
undefinedPoints = zeros(numProblemPoints, 3);

for count=1:numProblemPoints

    current = uniqueProblemNodes(count);
    
    rows = [find(bmesh(:,1) == current) find(bmesh(:,2) == current)];
    possibleNeighbours = [bmesh(rows,1) bmesh(rows,2)];    

    realNeighbours = (possibleNeighbours(possibleNeighbours ~= current));

    undefinedPoints(count,1) = current;
    undefinedPoints(count,2) = realNeighbours(1);
    undefinedPoints(count,3) = realNeighbours(2);

end


end


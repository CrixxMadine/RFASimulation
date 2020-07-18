function undefinedPoints = GetUndefinedBoundaryPoints(bmesh)

%% Function and parameter description

% Find problematic boundary nodes in grid
% 
% Problematic node happens if there are different boundary conditions
% colliding at one point e.g. the electrodes on the probe
%
% Since the dirichlet value has higher priority and therefore set the its value,
% the next neighbour on the neumann or robin side is calculated wrong
% 
% This function is to detect these problematic nodes
%
% Input arg:
% bmesh, where each row represents a boundary edge
% First and second column are the bnodes, third column is boundary type
% Boundary types: 1 -> Dirichlet, 2 -> Neumann, 3 -> Robin
%
% Returns:
% undefinedPoints: -> n x 3 matrix with the problematic node information
% First column contains the problematic point in pmesh
% Second and third column are the neigbour nodes of the undefined node


%% TODO clean up

rowsDirichCondition = find(bmesh(:,3) == 1);
allDirichNodes = unique([bmesh(rowsDirichCondition,1), bmesh(rowsDirichCondition,2)]);

rowsNeumannCondition = find(bmesh(:,3) ~= 1);
allNeumannNodes = unique([bmesh(rowsNeumannCondition,1), bmesh(rowsNeumannCondition,2)]);

% Node with both edges with dirichlet condition and other condition
doubleDefinedNodes = intersect(allDirichNodes, allNeumannNodes);

numDoubleDefinedPoints = size(doubleDefinedNodes,1);

% Case no problematic points (e.g. only one type of boundary conditions)
if (numDoubleDefinedPoints == 0) 
    undefinedPoints = zeros(0,0);
    return;
end 

problematicNodes(numDoubleDefinedPoints) = 0;

for count=1:length(doubleDefinedNodes)
corner = doubleDefinedNodes(count);
temp1 = find(bmesh(:,1)==corner);
temp2 = find(bmesh(:,2)==corner);
rowVonEcke = unique([temp1, temp2]);

bla = [bmesh(rowVonEcke,1) bmesh(rowVonEcke,2)];

yuhu = setdiff(bla, allDirichNodes);

% if (length(yuhu) > 0)
 problematicNodes(count) = yuhu;
%end

otherPoints = setdiff(bla,yuhu);
testPoints(count,1) = yuhu;
testPoints(count,2) = otherPoints(1);
testPoints(count,3) = otherPoints(2);

end

% On coarse grid it is possible, there is only one problem node 
% in between two dirichlit nodes, e.g. between the electrodes
% Given that case, number of problematic nodes reduces
temp = unique(problematicNodes);
uniqueProblemNodes = temp(temp~=0);

numProblemPoints = length(uniqueProblemNodes);
undefinedPoints = zeros(numProblemPoints, 3);

for i=1:numProblemPoints

current = uniqueProblemNodes(i);
% currentTest = problematicNodes(i);
rows = [find(bmesh(:,1) == current) find(bmesh(:,2) == current)];
neighbours = [bmesh(rows,1) bmesh(rows,2)];    
    
realNeighbours = (neighbours(neighbours ~= current));

undefinedPoints(i,1) = current;
undefinedPoints(i,2) = realNeighbours(1);
undefinedPoints(i,3) = realNeighbours(2);

end


end


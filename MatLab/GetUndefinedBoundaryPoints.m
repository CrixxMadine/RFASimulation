function undefinedPoints = GetUndefinedBoundaryPoints(bmesh)

hereIsDirich = find(bmesh(:,3) == 1);
theseAreDirichNodes = unique([bmesh(hereIsDirich,1), bmesh(hereIsDirich,2)]);

hereIsNeumann = find(bmesh(:,3) == 2);
neumannNodes = unique([bmesh(hereIsNeumann,1), bmesh(hereIsNeumann,2)]);

% Dirichlet corner nodes
[problemCornerNodes] = intersect(theseAreDirichNodes, neumannNodes);

numProblemPoints = size(problemCornerNodes,1);



% Case no problematic points (e.g. only one type of boundary conditions)
if (numProblemPoints == 0) 
    undefinedPoints = zeros(0,0);
    return;
end 

fuckYou(numProblemPoints) = 0;

for count=1:length(problemCornerNodes)
corner = problemCornerNodes(count);
temp1 = find(bmesh(:,1)==corner);
temp2 = find(bmesh(:,2)==corner);
rowVonEcke = unique([temp1, temp2]);

bla = [bmesh(rowVonEcke,1) bmesh(rowVonEcke,2)];

yuhu = setdiff(bla, theseAreDirichNodes);

if (length(yuhu) > 0)
fuckYou(count) = yuhu;
end

end

fast = unique(fuckYou);
fakeDir = fast(fast~=0);


undefinedPoints = zeros(numProblemPoints, 3);

for i=1:numProblemPoints

current = fakeDir(i);
rows = [find(bmesh(:,1) == current) find(bmesh(:,2) == current)];
neighbours = [bmesh(rows,1) bmesh(rows,2)];    
    
realNeighbours = (neighbours(neighbours ~= current));

undefinedPoints(i,1) = current;
undefinedPoints(i,2) = realNeighbours(1);
undefinedPoints(i,3) = realNeighbours(2);

end


end


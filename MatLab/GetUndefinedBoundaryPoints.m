function undefinedPoints = GetUndefinedBoundaryPoints(bmesh)

hereIsDirich = find(bmesh(:,3) == 1);
theseAreDirichNodes = unique([bmesh(hereIsDirich,1), bmesh(hereIsDirich,2)]);

hereIsNeumann = find(bmesh(:,3) == 2);
neumannNodes = unique([bmesh(hereIsNeumann,1), bmesh(hereIsNeumann,2)]);

% Dirichlet corner nodes
[vali] = intersect(theseAreDirichNodes, neumannNodes);

fuckYou(length(vali)) = 0;

for ff=1:length(vali)
ecke = vali(ff);
nnn = find(bmesh(:,1)==ecke);
uuu = find(bmesh(:,2)==ecke);
rowVonEcke = unique([nnn, uuu]);

bla = [bmesh(rowVonEcke,1) bmesh(rowVonEcke,2)];

yuhu = setdiff(bla, theseAreDirichNodes);

if (length(yuhu) > 0)
fuckYou(ff) = yuhu;
end

end

fast = unique(fuckYou);
fakeDir = fast(fast~=0);


numFakeDir = length(fakeDir);
undefinedPoints = zeros(numFakeDir,1));

for i=1:numFakeDir

current = fakeDir(i);
rows = [find(bmesh(:,1) == current) find(bmesh(:,2) == current)];
neighbours = [bmesh(rows,1) bmesh(rows,2)];    
    
realNeighbours = (neighbours(neighbours ~= current));

undefinedPoints(i,1) = current;
undefinedPoints(i,2) = realNeighbours(1);
undefinedPoints(i,3) = realNeighbours(2);

end


end


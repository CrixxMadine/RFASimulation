function undefinedPoints = GetUndefinedBoundaryPoints(bmesh)

hereIsDirich = find(bmesh(:,3) == 1);
theseAreDirichNodes = unique([bmesh(hereIsDirich,1), bmesh(hereIsDirich,2)]);

hereIsNeumann = find(bmesh(:,3) == 2);
neumannNodes = unique([bmesh(hereIsNeumann,1), bmesh(hereIsNeumann,2)]);

% Dirichlet corner nodes
[vali, posi] = intersect(theseAreDirichNodes, neumannNodes);

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


for i=1:length(fakeDir)


mmm = fakeDir(i);
rows = [find(bmesh(:,1) == mmm) find(bmesh(:,2) == mmm)];
neighbours = [bmesh(rows,1) bmesh(rows,2)];    
    
realNeighbours = (neighbours(neighbours ~= mmm));

uiuiu = realNeighbours(1);
jajaj = realNeighbours(2);

phi(mmm) = (phi(uiuiu) + phi(jajaj)) / 2;

end

undefinedPoints;

end


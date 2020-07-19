function [ux, uy] = VerticesGradientPDE(pmesh,tmesh, bedges, u)

% Calculate the gradient on the vertices from a PDE solution
% Gradient is undefined on boundaries by definition

ux = zeros(size(pmesh,1),1);
uy = zeros(size(pmesh,1),1);

[tx,ty] = pdegrad(pmesh',tmesh',u);


numTriangs = size(tmesh,1);
numNodes = size(pmesh,1);


allNodes = (1:numNodes);
boundNodes = unique([bedges(:,1), bedges(:,2)]);
innerNodes = setdiff(allNodes, boundNodes);

numInnerNodes = length(innerNodes);

for i=1:numInnerNodes
    node = innerNodes(i);
    
    % Find neighbors
    neighbourTriangs = find(tmesh == node);
    neighbourTriangs = mod(neighbourTriangs, numTriangs);
    
    % Possible error with modulo
    neighbourTriangs(neighbourTriangs==0) = numTriangs; 
    
    numNeighbours = length(neighbourTriangs);
    
    ux(node) = sum(tx(neighbourTriangs')) / numNeighbours;
    uy(node) = sum(tx(neighbourTriangs')) / numNeighbours;
          
end


tmesh_test = [tmesh, zeros(size(tmesh,1),1)];

figure(100);
pdeplot(pmesh',bedges',tmesh_test','xydata',tx,'zdata',ty)
title('Numerical gradient of the potential phi')
xlabel('r-axis');
ylabel('z-axis');

end


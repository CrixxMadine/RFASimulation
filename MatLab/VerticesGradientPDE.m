function [ux, uy] = VerticesGradientPDE(pmesh,tmesh, bedges, u)

% Calculate the gradient on the vertices from a PDE solution
% Gradient is undefined on boundaries by definition

ux = zeros(size(pmesh,1));
uy = zeros(size(pmesh,1));

[r,s] = pdegrad(pmesh',tmesh',u);

tx = r';
ty = s';


%allPoints = 


tmesh_test = [tmesh, zeros(size(tmesh,1),1)];

figure(300);
[tx,ty] = pdegrad(pmesh',tmesh',u);
pdeplot(pmesh',bedges',tmesh_test','xydata',tx,'zdata',ty)

end


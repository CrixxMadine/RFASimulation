function [sumOfElementIntegrals] = SurfaceIntegralTriangles(tmesh, pmesh, z_val)

% Return the surface integral on a triangulation
% 
% returns:
%  sum := numeric integral over the whole triangulation 
%
% Input args:
%  pmesh := point matrix of triangulation
%  tmesh := triangle matrix of triangulation
%  z_val := discrete values for every point of the triangulation


numberTriangles = size(tmesh,1);
sumOfElementIntegrals = 0;


for i=1:numberTriangles
   
    % get the global point numbers of the current triangle
    p = [tmesh(i,1),  tmesh(i,2), tmesh(i,3)];
    
    % get the coordinates of each point
    a = pmesh(p(1),:);
    b = pmesh(p(2),:);
    c = pmesh(p(3),:);
    
    % Area of triangle with given coordinates     
    A = [  1   ,  1   ,  1   ;
          a(1) , b(2) , c(2) ; 
          a(2) , b(2) , c(2) ];    
    
    area = 0.5 * abs(det(A));

    elementIntegral = 0.5 * (z_val(p(1)) + z_val(p(2)) + z_val(p(3))) * area;

    sumOfElementIntegrals = sumOfElementIntegrals + abs(elementIntegral);

end


end


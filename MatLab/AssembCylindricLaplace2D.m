function [Ah, fh] = AssembCylindricLaplace2D(pmesh, tmesh, k, q, f_rhs, intyp)

%% Function summary and arguments description

% Assemble the system of equations for cylindircal Laplace 2D
% Assemble each element, no boundary conditions yet
% Model equation: - div( k * grad(u) ) + q * u = f_rhs

% returns:
% Ah := FE-matrix (sum of stiffnessmatrix and q-massmatrix)
% fh := right hand side of the system of equations

% Input args:
% pmesh := point matrix of triangulation
% tmesh := triangle matrix of triangulation
% k, q  := coeffiients of the pde
% f_rhs := right hand side of the pde
%
% eltyp := type of regression function
% -> 1 = linear triangle element
%
% intyp := type of numerical quadrature
% -> 0 = with matlab library function
% -> 1 = own implementation of gauss quadrature


% explanations on the grid args (pmesh, emesh, tmesh):
%
% pmesh:  point matrix 
% Contains all the r and z coordinates of the points
%
% tmesh:  triangle matrix 
% First three lines: corners of the trinagle, counted counterclockwise 
% Fourth line: subdomain number (to be reserved;
%              not used by now, but could be used for material properties)
%
% bmesh:  boundray edge matrix : information on edges and its types 
%         -> to be used for boundary conditions


%% Define helper parameters

Ng = length(pmesh);    % total number of points
Me = length(tmesh);    % total number of triangle elements
Ne = 3;    % number of points per element - (3 for triangle)



%% Assemble the system of equations, elementwise

Ah = zeros(Ng,Ng);
fh = zeros(Ng,1);

for i=1:Me
    
    % get the global point numbers of the current triangle
    a = [tmesh(1,i),  tmesh(2,i), tmesh(3,i)];
    
    % get the coordinates of each point
    p1 = pmesh(:,a(1));
    p2 = pmesh(:,a(2));
    p3 = pmesh(:,a(3));
    
    % assemble matrix and right hand sand for the current element
    [K_elem, f_elem] = LinearElement2D(k, q, f_rhs, p1, p2, p3, intyp);
    
    % add element matrix and right hand side to the existing system
    for m = 1:Ne  
        
       fh(a(m)) = fh(a(m)) + f_elem(m);  
      
       for n = 1:Ne
           
           Ah(a(m),a(n)) = Ah(a(m),a(n)) + K_elem(m,n);
           
       end 
       
    end 
    
end 

end


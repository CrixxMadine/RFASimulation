function [Ah, fh] = AssembCylindricLaplace2D(pmesh, tmesh, k, q, f_rhs, intyp)

%% Function summary and arguments description

% Assemble the system of equations for cylindircal Laplace 2D
% Assemble each element, no boundary conditions yet
% Model equation: - div( k * grad(u) ) + q * u = f_rhs

% returns:
%  Ah := sparse FE-matrix (sum of stiffnessmatrix and q-massmatrix)
%  fh := right hand side of the system of equations

% Input args:
%  pmesh := point matrix of triangulation
%  tmesh := triangle matrix of triangulation
%  k, q  := coeffiients of the pde
%  f_rhs := right hand side of the pde
%
%  eltyp := type of regression function
%  -> 1 = linear triangle element
%
%  intyp := type of numerical quadrature
%  -> 0 = with matlab library function
%  -> 1 = own implementation of gauss quadrature


%  explanations on the grid args (pmesh, emesh, tmesh):
%
%  pmesh:  point matrix 
%  Contains all the r and z coordinates of the points
%
%  tmesh:  triangle matrix 
%  First three lines: corners of the trinagle, counted counterclockwise 
%  Fourth line: subdomain number (to be reserved;
%               not used by now, but could be used for material properties)
%
%  bmesh:  boundray edge matrix : information on edges and its types 
%         -> to be used for boundary conditions


%% Implementation


% Define helper parameters

Ng = length(pmesh);    % total number of points
Me = length(tmesh);    % total number of triangle elements
Ne = 3;    % number of points per element - (3 for triangle)


% Assemble the system of equations, elementwise

Ah_row = zeros(3*Ng,1);
Ah_col = Ah_row;
Ah_val = Ah_row;

fh = zeros(Ng,1);
count = 0;

for i=1:Me
    
    % Get the global point numbers of the current triangle
    a = [tmesh(i,1),  tmesh(i,2), tmesh(i,3)];
    
    % Get the coordinates of each point
    p1 = pmesh(a(1),:)';
    p2 = pmesh(a(2),:)';
    p3 = pmesh(a(3),:)';
    
    % Assemble matrix and right hand sand for the current element
    [K_elem, f_elem] = LinearElement2D(k, q, f_rhs, p1, p2, p3, intyp);
    
    % Add element matrix and right hand side to the existing system
    for m = 1:Ne  
        
       fh(a(m)) = fh(a(m)) + f_elem(m);  
      
       for n = 1:Ne
                     
           count = count + 1;
                      
           Ah_row(count) = a(m);
           Ah_col(count) = a(n);
           Ah_val(count) = K_elem(m,n); 
           
       end 
       
    end 
    
end 

Ah = sparse(Ah_row, Ah_col, Ah_val);

end


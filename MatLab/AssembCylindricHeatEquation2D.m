function [Ah, Mh, fh] = AssembCylindricHeatEquation2D(pmesh, tmesh, k, q, Q_val, intyp)

%% Function summary and arguments description

% Assemble the system of equations for cylindircal Laplace 2D
% Model equation: du/dt - div( k * grad(u) ) + q * u = f_rhs

% Assemble each element, no boundary conditions yet

% returns:
% Ah := FE-matrix (sum of stiffnessmatrix and q-massmatrix)
% Mh := M-Matrix (the prefactor matrix of du/dt)
% fh := right hand side of the system of equations

% -- Inout args are very similar to AssembleCylindricLaPlace --
% New arg: Q_val -> replacement for f_rhs
% Discrete vlaue for the right hand side

%% Define additional parameters

Ng = length(pmesh);    % total number of points
Me = length(tmesh);    % total number of triangle elements
Ne = 3;    % number of points per element - (3 for triangle)


%% Assemble the system of equations, elementwise

Ah = zeros(Ng,Ng);
Mh = zeros(Ng,Ng);
fh = zeros(Ng,1);

for i=1:Me
    
    % get the global point numbers of the current triangle
    a = [tmesh(1,i),  tmesh(2,i), tmesh(3,i)];
    
    % get the coordinates of each point
    p1 = pmesh(:,a(1));
    p2 = pmesh(:,a(2));
    p3 = pmesh(:,a(3));
    
    % Define f_rhs as constant average value from the triangulation
    f_rhs = @(r,z) 1 / 3 * (Q_val(a(1)) + Q_val(a(2)) + Q_val(a(3)));
    
    % assemble matrix and right hand sand for the current element
    [K_elem, f_elem, M_elem] = LinearElement2D(k, q, f_rhs, p1, p2, p3, intyp, 'true');
    
    % add element matrix and right hand side to the existing system
    for m = 1:Ne  
        
       fh(a(m)) = fh(a(m)) + f_elem(m);  
      
       for n = 1:Ne
           
           Ah(a(m),a(n)) = Ah(a(m),a(n)) + K_elem(m,n);
           
           Mh(a(m),a(n)) = Mh(a(m),a(n)) + M_elem(m,n);
       end 
       
    end 
    
end 

end



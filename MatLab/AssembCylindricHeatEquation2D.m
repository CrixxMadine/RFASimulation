function [Ah, Mh, fh] = AssembCylindricHeatEquation2D(pmesh, tmesh, k, q, Q_rfa, Q_perf, intyp)

%% Function summary and arguments description

% Assemble the system of equations for cylindircal Laplace 2D
% Model equation: du/dt - div( k * grad(u) ) + q * u = f_rhs

% Assemble each element, no boundary conditions yet

% returns:
% Ah := sparse FE-matrix (sum of stiffnessmatrix and q-massmatrix)
% Mh := sparse M-Matrix (the prefactor matrix of du/dt)
% fh := right hand side of the system of equations

% -- Inout args are very similar to AssembleCylindricLaPlace --
% New arg: Q_val -> replacement for f_rhs
% Discrete vlaue for the right hand side

%% Define additional parameters

Ng = length(pmesh);    % total number of points
Me = length(tmesh);    % total number of triangle elements
Ne = 3;    % number of points per element - (3 for triangle)


%% Assemble the system of equations, elementwise

Ah_row = zeros(3*Ng,1);
Ah_col = Ah_row;
Ah_val = Ah_row;

Mh_row = zeros(3*Ng,1);
Mh_col = Mh_row;
Mh_val = Mh_row;

fh = zeros(Ng,1);

count = 0;

% OLD -> Using full matrix
% Ah = zeros(Ng,Ng);
% Mh = zeros(Ng,Ng);

for i=1:Me
    
    % get the global point numbers of the current triangle
    a = [tmesh(i,1),  tmesh(i,2), tmesh(i,3)];
    
    % get the coordinates of each point
    p1 = pmesh(a(1),:)';
    p2 = pmesh(a(2),:)';
    p3 = pmesh(a(3),:)';
    
    % Define f_rhs as constant average value from the triangulation
    % f_rhs = @(r,z) 1 / 3 * (Q_val(a(1)) + Q_val(a(2)) + Q_val(a(3)));
    
    Q_perf_val = 1 / 3 * Q_perf(a(1)) + Q_perf(a(2)) + Q_perf(a(3));
    % Better use pdegrad for triangle
    f_rhs = @(r,z) Q_rfa(i) + Q_perf_val;
    
    % assemble matrix and right hand sand for the current element
    [K_elem, f_elem, M_elem] = LinearElement2D(k, q, f_rhs, p1, p2, p3, intyp, 'true');
    
    % add element matrix and right hand side to the existing system
    for m = 1:Ne  
        
       fh(a(m)) = fh(a(m)) + f_elem(m);  
      
       for n = 1:Ne

           count = count + 1;
                      
           Ah_row(count) = a(m);
           Ah_col(count) = a(n);
           Ah_val(count) = K_elem(m,n); 
           
           Mh_row(count) = a(m);
           Mh_col(count) = a(n);
           Mh_val(count) = M_elem(m,n); 
           
           % OLD -> Using full matrix
           % Ah(a(m),a(n)) = Ah(a(m),a(n)) + K_elem(m,n);          
           % Mh(a(m),a(n)) = Mh(a(m),a(n)) + M_elem(m,n);
           
       end 
       
    end 
    
end 

Ah = sparse(Ah_row, Ah_col, Ah_val);
Mh = sparse(Mh_row, Mh_col, Mh_val);

end



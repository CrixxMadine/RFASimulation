function [K, f, M] = LinearElement2D(k, q, f_rhs, a1, a2, a3, intyp, opt_M)

%% Function and input argument discription

% Calculates a single linear 2D triangle element with linear basis function
% Uses cylinder coordinates for calculation of integrals
% This function is to be seen as a documentation of the algorithm
% It is not optimized for complexity and runtime
% For an optimized simulation see C++ Code

% Returns equation : M * du/dt - K * u = f

% Returns
% K := block matrix of the element (stiffness matrix + q_matrix)
% f := right side of the element
% M := block parabolic mass matrix of the element (optional)

% Iput arguments
% k, q  := coefficients of PDE    ->  arg is vector
% f_rhs := right side of the PDE  ->  arg is vector
%
% a1 = (r, z) := coordinates of the left corner of the triangle
% a2,a3       := other corners counted counterclockwise
%
% intyp = integration type
% for reference see 'QuadratureTriangle2D'
% 
% opt_M := optionally return prefactor matrix of du/dt in parabolic 

if nargin==7      % you must specify if you want M for parabolic returned
   opt_M = 'false';
elseif nargin<7
	error('You did not provide enough input arguments.')
elseif nargin == 8
    if ~ischar(opt_M)
        error('Last input argument must be a string with value true or false')
    end
end



%% Transformation on general reference triangle [(0,0), (0,1), (1,0)]

% for ref see script 'Steinbach Numerik 3', chap.: 3.4.2.1 and 3.4.3

J = [(a2-a1), (a3-a1)];  % jacobi matrix
                      
JT_inv = inv(J)';        % transposed of the inverted jacobi matrix
                  
det_J = abs(det(J));     % absolute value of determinant
                  

r_ref = [JT_inv(1,1) JT_inv(1,2)];  % partial jacobi matrix dr
z_ref = [JT_inv(2,1) JT_inv(2,2)];  % partial jacobi matrix dz


if (intyp == 1 || intyp == 2)
    
    F = @(my,ny) mtimes(J, [my; ny]) + a1;  
  
    
end

% these are the corners of normalized reference triangular
a = [0; 0];
b = [1; 0];
c = [0; 1];


%% Define linear basis functions

phi1 = @(my, ny)  1 - my - ny;   % phi1
phi2 = @(my, ny)  my;            % phi2
phi3 = @(my, ny)  ny;            % phi3

dr_phi1 = @(my,ny) -1;
dr_phi2 = @(my,ny)  1;
dr_phi3 = @(my,ny)  0;

dz_phi1 = @(my,ny) -1;
dz_phi2 = @(my,ny)  0;
dz_phi3 = @(my,ny)  1;


%% Define integrals for quadrature and associated right hand side

%% Cylindrical LaPlace Equation

% Equation: \integrate: (du/dr * dv/dr + du/dz * dv/dz) * r * dr dz
cyl_int_11 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi1(r,z) , r_ref .* dr_phi1(r,z)) + ... 
                               dot(z_ref .* dz_phi1(r,z) , z_ref .* dz_phi1(r,z))) * r;

cyl_int_12 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi1(r,z) , r_ref .* dr_phi2(r,z)) + ... 
                               dot(z_ref .* dz_phi1(r,z) , z_ref .* dz_phi2(r,z))) * r;
                           
cyl_int_13 = @(r,z) det_J  * k(F(r,z)) .* ( dot(r_ref .* dr_phi1(r,z) , r_ref .* dr_phi3(r,z)) + ... 
                               dot(z_ref .* dz_phi1(r,z) , z_ref .* dz_phi3(r,z))) * r;                           
                           
cyl_int_21 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi2(r,z) , r_ref .* dr_phi1(r,z)) + ... 
                               dot(z_ref .* dz_phi2(r,z) , z_ref .* dz_phi1(r,z))) * r;

cyl_int_22 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi2(r,z) , r_ref .* dr_phi2(r,z)) + ... 
                               dot(z_ref .* dz_phi2(r,z) , z_ref .* dz_phi2(r,z))) * r;
                           
cyl_int_23 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi2(r,z) , r_ref .* dr_phi3(r,z)) + ... 
                               dot(z_ref .* dz_phi2(r,z) , z_ref .* dz_phi3(r,z))) * r;                                
                           
cyl_int_31 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi3(r,z) , r_ref .* dr_phi1(r,z)) + ... 
                               dot(z_ref .* dz_phi3(r,z) , z_ref .* dz_phi1(r,z))) * r;

cyl_int_32 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi3(r,z) , r_ref .* dr_phi2(r,z)) + ... 
                               dot(z_ref .* dz_phi3(r,z) , z_ref .* dz_phi2(r,z))) * r;
                           
cyl_int_33 = @(r,z) det_J .* k(F(r,z)) .* ( dot(r_ref .* dr_phi3(r,z) , r_ref .* dr_phi3(r,z)) + ... 
                               dot(z_ref .* dz_phi3(r,z) , z_ref .* dz_phi3(r,z))) * r;
                             
cyl_int = {cyl_int_11, cyl_int_12, cyl_int_13;
           cyl_int_21, cyl_int_22, cyl_int_23;
           cyl_int_31, cyl_int_32, cyl_int_33};

% Equation: \integrate: f * v * r * dr dz    
cyl_fk_int1 = @(r,z) det_J .* f_rhs(F(r,z)) .* phi1(r,z) .* r; 
cyl_fk_int2 = @(r,z) det_J .* f_rhs(F(r,z)) .* phi2(r,z) .* r;
cyl_fk_int3 = @(r,z) det_J .* f_rhs(F(r,z)) .* phi3(r,z) .* r;

cyl_fk_int = {cyl_fk_int1 ; cyl_fk_int2 ; cyl_fk_int3};


%% Optional mass matrix for parabolic equations

if strcmp(opt_M, 'true')
    
    % Equation: \integrate: u * v * r * dr dz
    mass_int_11 = @(r,z) det_J .* phi1(r,z) .* phi1(r,z);
    mass_int_12 = @(r,z) det_J .* phi1(r,z) .* phi2(r,z);
    mass_int_13 = @(r,z) det_J .* phi1(r,z) .* phi3(r,z);

    mass_int_21 = @(r,z) det_J .* phi2(r,z) .* phi1(r,z);
    mass_int_22 = @(r,z) det_J .* phi2(r,z) .* phi2(r,z);
    mass_int_23 = @(r,z) det_J .* phi2(r,z) .* phi3(r,z);

    mass_int_31 = @(r,z) det_J .* phi3(r,z) .* phi1(r,z);
    mass_int_32 = @(r,z) det_J .* phi3(r,z) .* phi2(r,z);
    mass_int_33 = @(r,z) det_J .* phi3(r,z) .* phi3(r,z);

    mass_int = {mass_int_11, mass_int_12, mass_int_13;
                mass_int_21, mass_int_22, mass_int_23;
                mass_int_31, mass_int_32, mass_int_33};
end        


%% Calculate the element matrix and associated right hand side

K = zeros(3,3);
M = zeros(3,3);
f = zeros(3,1);

for alpha=1:3
    
    for beta=1:3      
        K(alpha, beta) = QuadratureTriangle2D(cyl_int{alpha, beta}, intyp, a, b, c);
        
        if strcmp(opt_M, 'true')
            M(alpha, beta) = QuadratureTriangle2D(mass_int{alpha, beta}, intyp, a, b, c);
        end
    end
    
    f(alpha) = QuadratureTriangle2D(cyl_fk_int{alpha, 1}, intyp, a, b, c);
    
end



end


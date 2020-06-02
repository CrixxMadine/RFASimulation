function [K, f] = LinearElement2D(k, q, f_rhs, a1, a2, a3, intyp)

%% Function and input argument discription

% Calculates a single linear 2D triangle element with linear basis function
% This function is to be seen as a documentation of the algorithm
% It is not optimized for complexity and runtime
% For an optimized simulation see C++ Code

% Returns
% K := block matrix of the element
% f := right side of the element

% Iput arguments
% k, q  := coefficients of PDE    ->  arg is vector
% f_rhs := right side of the PDE  ->  arg is vector
%
% a1 = (r, z) := coordinates of the left corner of the triangle
% a2,a3       := other corners counted counterclockwise
%
% intyp = integration type
% for reference see 'QuadratureTriangle2D'


%% Transformation on general reference triangle [(0,0), (0,1), (1,0)]

% for ref see script 'Steinbach Numerik 3', chap.: 3.4.2.1 and 3.4.3

J = [(a2-a1), (a3-a1)];  % jacobi matrix
                      
JT_inv = inv(J)';        % transposed of the inverted jacobi matrix
                  
det_J = abs(det(J));     % absolute value of determinant
                  

if (intyp == 1 || intyp == 2)
    
    F = @(my,ny) mtimes(J, [my; ny]) + a1;  

else
    
    % matlab library integration uses scalar operations
    % function F_scalar is defined at the end of this script
    F = @F_scalar;    
    
end

% these are the corners of normalized reference triangular
a = [0; 0];
b = [1; 0];
c = [0; 1];


%% Define linear basis functions

phi1 = @(my, ny)  1 - my - ny;   % phi1
phi2 = @(my, ny)  my;            % phi2
phi3 = @(my, ny)  ny;            % phi3

r_ref = 1 / (a2(1) - a1(1));
z_ref = 1 / (a3(2) - a1(2));

dr_phi1 = @(my,ny) -1;
dr_phi2 = @(my,ny)  1;
dr_phi3 = @(my,ny)  0;

dz_phi1 = @(my,ny) -1;
dz_phi2 = @(my,ny)  0;
dz_phi3 = @(my,ny)  1;

% note: gradient is defined as grad f(x,y) = [df/dx, df/dy]
grad_phi1 = @(my, ny)  [-1 -1];  % grad phi1
grad_phi2 = @(my, ny)  [ 1  0];  % grad phi2 
grad_phi3 = @(my, ny)  [ 0  1];  % grad phi3 


%% Define integrals for quadrature and associated right hand side

% remember, this is an academic approach
% optimized values are implemented in C++ algorithm

%% Cylindrical LaPlace Equation

% Equation: \integrate: (dr u * dr v + dz u * dz v ) * r * dr dz

cyl_int_11 = @(r,z) det_J .* (r_ref * dr_phi1(r,z) * r_ref * dr_phi1(r,z) + ... 
                              z_ref * dz_phi1(r,z) * z_ref * dr_phi1(r,z)) * r;
cyl_int_12 = @(r,z) det_J .* (r_ref * dr_phi1(r,z) * r_ref * dr_phi2(r,z) + ... 
                              z_ref * dz_phi1(r,z) * z_ref * dr_phi2(r,z)) * r; 
cyl_int_13 = @(r,z) det_J .* (r_ref * dr_phi1(r,z) * r_ref * dr_phi3(r,z) + ... 
                              z_ref * dz_phi1(r,z) * z_ref * dr_phi3(r,z)) * r; 

cyl_int_21 = @(r,z) det_J .* (r_ref * dr_phi2(r,z) * r_ref * dr_phi1(r,z) + ... 
                              z_ref * dz_phi2(r,z) * z_ref * dr_phi1(r,z)) * r;
cyl_int_22 = @(r,z) det_J .* (r_ref * dr_phi2(r,z) * r_ref * dr_phi2(r,z) + ... 
                              z_ref * dz_phi2(r,z) * z_ref * dr_phi2(r,z)) * r; 
cyl_int_23 = @(r,z) det_J .* (r_ref * dr_phi2(r,z) * r_ref * dr_phi3(r,z) + ... 
                              z_ref * dz_phi2(r,z) * z_ref * dr_phi3(r,z)) * r; 

cyl_int_31 = @(r,z) det_J .* (r_ref * dr_phi3(r,z) * r_ref * dr_phi1(r,z) + ... 
                              z_ref * dz_phi3(r,z) * z_ref * dr_phi1(r,z)) * r;
cyl_int_32 = @(r,z) det_J .* (r_ref * dr_phi3(r,z) * r_ref * dr_phi2(r,z) + ... 
                              z_ref * dz_phi3(r,z) * z_ref * dr_phi2(r,z)) * r; 
cyl_int_33 = @(r,z) det_J .* (r_ref * dr_phi3(r,z) * r_ref * dr_phi3(r,z) + ... 
                              z_ref * dz_phi3(r,z) * z_ref * dr_phi3(r,z)) * r; 
                   
cyl_int = {cyl_int_11, cyl_int_12, cyl_int_13;
           cyl_int_21, cyl_int_22, cyl_int_23;
           cyl_int_31, cyl_int_32, cyl_int_33};
       
       
%% 2D FEM model problem from Numerik 3, Exercise 3

integral_11 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi1(my,ny)') , (JT_inv * grad_phi1(my,ny)')) ...
                     + q(F(my,ny)) .* phi1(my,ny) .* phi1(my,ny) ) );
                                                                                  
integral_12 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi2(my,ny)') , (JT_inv * grad_phi1(my,ny)')) ...
                     + q(F(my,ny)) .* phi2(my,ny) .* phi1(my,ny) ) );

integral_13 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi3(my,ny)') , (JT_inv * grad_phi1(my,ny)')) ...
                     + q(F(my,ny)) .* phi3(my,ny) .* phi1(my,ny) ) );


integral_21 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi1(my,ny)') , (JT_inv * grad_phi2(my,ny)')) ...
                     + q(F(my,ny)) .* phi1(my,ny) .* phi2(my,ny) ) );

integral_22 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi2(my,ny)') , (JT_inv * grad_phi2(my,ny)')) ...
                     + q(F(my,ny)) .* phi2(my,ny) .* phi2(my,ny) ) );

integral_23 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi3(my,ny)') , (JT_inv * grad_phi2(my,ny)')) ...
                     + q(F(my,ny)) .* phi3(my,ny) .* phi2(my,ny) ) );
                 

integral_31 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi1(my,ny)') , (JT_inv * grad_phi3(my,ny)')) ...
                     + q(F(my,ny)) .* phi1(my,ny) .* phi3(my,ny) ) );

integral_32 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi2(my,ny)') , (JT_inv * grad_phi3(my,ny)')) ...
                     + q(F(my,ny)) .* phi2(my,ny) .* phi3(my,ny) ) );

integral_33 = @(my,ny) ( det_J .* ( k(F(my,ny)) .* dot((JT_inv * grad_phi3(my,ny)') , (JT_inv * grad_phi3(my,ny)')) ...
                     + q(F(my,ny)) .* phi3(my,ny) .* phi3(my,ny) ) );             


integrals = {integral_11, integral_12, integral_13;
             integral_21, integral_22, integral_23;
             integral_31, integral_32, integral_33};

         
fk_int1 = @(my,ny) det_J .* f_rhs(F(my,ny)) .* phi1(my,ny); 
fk_int2 = @(my,ny) det_J .* f_rhs(F(my,ny)) .* phi2(my,ny);
fk_int3 = @(my,ny) det_J .* f_rhs(F(my,ny)) .* phi3(my,ny);

fk_int = {fk_int1 ; fk_int2 ; fk_int3};


%% Calculate the element matrix and associated right hand side

K = zeros(3,3);
f = zeros(3,1);

for alpha=1:3
    
    for beta=1:3
        % Old Version
        % K(alpha, beta) = QuadratureTriangle2D(integrals{alpha, beta}, intyp, a, b, c);
        
        % Cyl LaPlace
        K(alpha, beta) = QuadratureTriangle2D(cyl_int{alpha, beta}, inty, a, b, c);
    end
    
    % Old Version
    % f(alpha) = QuadratureTriangle2D(fk_int{alpha, 1}, intyp, a, b, c);
    
    % Cyl LaPlace
    % Value is zero, so nothing o do here
    
end



function res = F_scalar(my,ny)

    % helper function to make F(x) scalar
    % mandatory for using matlab library function integral12 
    
    x = J(1,1) .* my + J(1,2) .* ny + a1(1);
    
    y = J(2,1) .* my + J(2,2) .* ny + a1(2);
    
    res = [x ; y];
    
end

end


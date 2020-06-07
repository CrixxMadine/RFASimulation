% This Skript is for experimenting with the RFA DGLs

% for reference see article of my advisor:
% Tim Kroeger et al. 
% "Numerical Simulation of Radio Frequency
%  Ablation with State Dependent Material
%  Parameters in Three Space Dimensions"

clear variables;
disp('MatLab ..... yeah :(')

%% Space and time dimensions

% cylindric domain  
%  ___   ___
% |   | |   |
% |   |_|   |
% |         |
% |_________|
%

% Transformation to zylindric coordinates 
x = @(r,phi,z) r * cos(phi);
y = @(r,phi,z) r * sin(phi);
z = @(r,phi,z) z;

% time dimension
t_start = 0.00;
t_step  = 0.01;
t_end   = 5.00;   
t_vec   = t_start:t_step:t_end;


%% Material parameters
% For the first run, everything should be absolutely const.

% Constant material parameters 
sigma  = 0; % electric conductivity
rho    = 0; % density
c      = 0; % specific heat capacity
lambda = 0; % thermal conductivty

% TODO

% Alternative model for linear dependent material parameters
% for ref. see books from Stein T. and Kröger T.
sigma  = [ 0.21     0.013      -1     1.143 ]'; 
rho    = [ 1080   -0.00056   -0.657     0   ]';
c      = [ 3455       0      -0.376     0   ]'; 
lambda = [ 0.437    0.0025      0       0   ]'; 

% variable state parameters
phi   = [ ]; % electric potential  
temp  = [ ]; % temperature distribution
F_wat = [ ]; % relative content of fluid water
F_vap = [ ]; % relative content of vapor
F_coa = [ ]; % coagulation state

% TODO 
% F_wat(t,x) + F_vap(t,x) == 1;  -> this would be intuitive
% F_wat(t,x) + F_vap(t,x) <= 1;  -> this is more appropriate

% Constant Prefactors for PDE
nu = @(x) 0.01765;  % Constant in first try


%% TRYING REGION SOLVE ELECTRIC POTENTIAL

% Grid Generation 

% Extra coarse grid 
% [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraCoarse\');
%  bmesh = DefineBoundaryConditions(bedges);
 
% Extra fine grid
[pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraFine\');
 bmesh = DefineBoundaryConditions(bedges);

% Plot the mesh, for control
figure(1);
subplot(1,2,1);
trimesh(tmesh', pmesh(1,:)', pmesh(2,:)');
title('The used triangulation of the domain');

subplot(1,2,2);
scatter(pmesh(1,:), pmesh(2,:));
title('All the points in the triangulation');


% Assemble PDE for electric potential
[Ah, fh] = AssembCylindricLaplace2D(pmesh, tmesh);

% Add boundary conditions
[Ah, fh] = AddBoundaryConditionsToFEMatrix(Ah, fh, pmesh, tmesh, bmesh);

% Solve the system of equations
uh = Ah \ fh;

% Calculate the power from the electric potential

% Calclulate power(x,y) ...
power = zeros(size(uh,1),1);

[phi_dx, phi_dy] = TriangularGradient(tmesh, pmesh, uh);

for i=1:size(power,1)
    power(i) = sigma(1) * norm([phi_dx(i), phi_dy(i)])^2;   
end

% Calculate total power
totalPower = SurfaceIntegralTriangles(tmesh, pmesh, power);

thisIsJustForBrakePoint = 0;

figure(2);
trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', power);
title('Solution of the finite element method');


stopTheExecutionHereBreakpoint = 0;



%% Information on solving the PDE's

% In Numerik 3 an elliptical PDE problem was discussed

% We can solve the classic problem as 
% -div( k(x)) grad u(x) ) + q(x) u(x) = f (x)

% x can be a vector of any dimension up to 3 dimensions
% every boundary condition is allowed, with a constraint on Robin

% With the following constraints:
% assuming k(x) ELEMENT OF C^1{OMEGA}
% assuming q(x), f(x) ELEMENT OF C{OEMGA}
% assuming k(x) >= k_0 > 0  (must be positive)
% assuming q(x) >= 0        (can be positive or zero) 
% -> In case of robin boundary conditions: 
% assuming prefactor kappa(x) >= 0 


%% Electric energy (EE)

% TODO

% We can use the elliptical approach above to compute phi(t,x)
% Assuming phi(t,x) is quasistatic, we can model it as follows

%% phi: Fixed potential on electrodes (I x OMEGA_elec)

% -- First try only, this is very arbitrarily
% -- Actually the potential is induced by a function
phi_electrode1 = +1;
phi_electrode2 = -1;


%% phi: Inner domain (I x (OMEGA \ OMEGA_elec) )

% In the inner domain, phi is modeled as follows
% - div ( (sigma(x,y,z) * grad phi(x,y,z) ) = 0

% this is an elliptic problem analogous to Numerik 3
% q(x) = 0, so the mass matrix is empty

% If sigma is constant, you can rewrite the problem
% - sigma * LAPLACE ( phi(x,y,z) ) = 0
% So it becomes type of Laplace's Equation

% Equation in cylindric coordinates
% - 1/r dphi/dr - d^2phi/r^2 - d^2phi/dz^2 = 0 

k_inner_domain = 0;  % Assuming sigma is constant
q_inner_domain = 0;
f_inner_domain = 0;

% TODO: this is phi, calculate it
u_inner_domain = [];   


%% phi: Robin boundary condition (I x GAMMA)
% n * NABLA phi = ( (n * (s-x)) / ( |s-x|^2 ) ) * phi
% s = barycenter of the union of all probes

% TODO

%% E: calculate the power

% TODO

% In case of constant material parameters:
% Phi becomes linear and time independent
% Since phi is time independent, we only have to solve it once

% Q_rf(t,x) = p((t,x)) * (p_eff(t)/p_total(t))



%% Temperature  
Q_rf   = @(t,x) 1;     % 
% TODO

Q_perf = @(t,x) 1;     % Q_perfusion (~Durchströmung)
% TODO    

Q_pc   = @(t,x) 0;     % Q_phase_change
% Is not considered in first simulation 
% Cause its model not reliably (Kroeger)

pde_phi = @(x) 1;
pde_Q = @(t,x) Q_rf(t,x) + Q_perf(t,x) + Q_pc(t,x);



%% Abbreviations
%
% ----- Keywords ------
% Keywords in comments are written in capital letters
% 
% NABLA  : The nabla operator
%
% DIV    : divergence,       NABLA /dotproduct something
% GRAD   : gradient,         NABLA /multiply something
% LAPLACE: Laplace operator, NABLA /dotproduct NABLA /multiply something
%
% PARAM  : this marks a material parameter
% PDE    : this marks a PDE
% STATE  : this marks a depending state
% 
%


% ----- All Abbreviations in order -----
%
% c      : PARAM - specific heat capacity
%  
% F      : marks stated related to fluid
% F_wat  : STATE - relative content of fluid water
% F_vap  : STATE - relative content of vapor
% F_coa  : STARE - coagulation state
%
% gamma  : the outer boundary of the domain
%
% lambda : PARAM - thermal conductivity
% 
% omega  : the whole domain
%
% phi    : PDE, STATE - electric potential  phi(t,x)
%
% Q      : PDE  -  heat energy Q(t,x)
% 
% Q_pc   : phase change part of Q
% Q_perf : perfusion part of Q 
% Q_rf   : radio frequency part of Q
%
% rho    : PARAM - density
%
% sigma  : PARAM - electric conductivity
%
% Temp   : STATE - temperature distribution
%
%
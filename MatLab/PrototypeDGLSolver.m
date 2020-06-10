
%% General information about the script

% The whole matlab code is part of an AMP application project
% Project is about numeric simulation of radio frequency ablation

% This skript is mainly for experimenting with the RFA PDEs
% Also general documentation of the numeric approaches

% THIS SCRIPT IS BY NO MEANS OPTIMIZED
% It massively prefers intention telling over optimization
% For an optimized simulation program, see the C++ source code

% for model reference see article of my advisor:
% Tim Kroeger et al. 
% "Numerical Simulation of Radio Frequency
%  Ablation with State Dependent Material
%  Parameters in Three Space Dimensions"

%  TODO: ADD PUBLISHER

clear variables;
disp('MatLab ..... yeah :(')

%% Space and time dimensions

% This simulation models a needle inserted in a malignant tissue
% Due to axis symmetrie, the whole simulation is reduced to a 2D-problem
% All calculations are done in cylindrical coordinates

% cylindrical domain for the simulation
%  ____    ____
% |    |  |    |
% |    |  |    |
% |    |  |    |
% |     \/     |
% |            |
% |            |
% |____________|
%

% Transformation reference to cylindrical coordinates 
x = @(r,phi,z) r * cos(phi);
y = @(r,phi,z) r * sin(phi);
z = @(r,phi,z) z;

% Time discretization in seconds for time-dependant simulation 
t_start = 0.00;  % Starting point -> t = 0 seconds
t_step  = 0.01;
t_end   = 5.00;  
t_vec   = t_start:t_step:t_end;


%% Material parameters

% For the first simulation run, everything should be absolutely const.

% Alternative model for linear dependent material parameters
% for ref. see books from Stein T. and Kröger T.
sigma_stein  = [ 0.21     0.013      -1     1.143 ]'; 
rho_stein    = [ 1080   -0.00056   -0.657     0   ]';
c_stein      = [ 3455       0      -0.376     0   ]'; 
lambda_stein = [ 0.437    0.0025      0       0   ]'; 

% Constant material parameters 
sigma_phi = sigma_stein(1);  % electric conductivity
rho_blood = rho_stein(1);    % density
c_blood   = c_stein(1);      % specific heat capacity
lambda    = lambda_stein(1); % thermal conductivty

nu_blood  = 0.01765;         % blood perfusion coefficient

% Variable state parameters
phi   = [ ]; % electric potential  
temp  = [ ]; % temperature distribution
F_wat = [ ]; % relative content of fluid water
F_vap = [ ]; % relative content of vapor
F_coa = [ ]; % coagulation state

% TODO 
% F_wat(t,x) + F_vap(t,x) == 1;  -> this would be intuitive
% F_wat(t,x) + F_vap(t,x) <= 1;  -> this is more appropriate

% Constant Prefactors for PDE



%% Grid Generation 

% Extra coarse grid 
% [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraCoarse\');
%  bmesh = DefineBoundaryConditions(bedges);
 
% Extra fine grid
[pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraFine\');

% Plot the mesh, for control
figure(1);
subplot(1,2,1);
trimesh(tmesh', pmesh(1,:)', pmesh(2,:)');
title('The used triangulation of the domain');

subplot(1,2,2);
scatter(pmesh(1,:), pmesh(2,:));
title('All the points in the triangulation');

%% Calculate electrical potential phi 

% In the inner domain, phi is quasistatic and modeled as follows
% |   - div ( (sigma(x,y,z) * grad phi(x,y,z) ) = 0   |
%  -> This is a Laplacian Equation, elliptical PDE second order

% Add boundary conditions for the elliptical problem
bmesh = DefineBoundaryConditions(bedges, 'phi');

% Define specific parameters for phi PDE

k_EPot  = @(r,z) 1;
q_EPot  = @(r,z) 0;
f_rhs_EPot = @(r,z) 0;
intyp = 1;

% Assemble FEM system of equations 
[Ah, fh] = AssembCylindricLaplace2D(pmesh, tmesh, k_EPot, q_EPot, f_rhs_EPot, intyp);

% Add boundary conditions
[Ah, fh] = AddBoundaryConditionsToFEMatrix(Ah, fh, pmesh, tmesh, bmesh);

% Solve the system of equations
uh = Ah \ fh;

figure(2);
trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', uh);
title('Solution of the finite element method for phi');

%% Calculate electric power from the electric potential

power = zeros(size(uh,1),1);

% Get the numerical gradient of every vertex
[phi_dx, phi_dy] = TriangularGradient(tmesh, pmesh, uh);

% Calclulate power(r,z) for every vertex

for i=1:size(power,1)
    power(i) = sigma_phi * norm([phi_dx(i), phi_dy(i)])^2;   
end

% Calculate total power of the domain
totalPower = SurfaceIntegralTriangles(tmesh, pmesh, power);

% Calculate effective power of the model 
power_setup = 200;   % power of the generator (in range 20-200 W)
U_elec = 2;          % Potential difference of the two electrodes

R_setup = 1; % TODO  % inner resistance of the generator
R_tis = U_elec * U_elec / totalPower; % tissue resistance

effectivePower = (4 * power_setup * R_tis * R_setup) / ...
    ((R_tis + R_setup)^2);  % effective power of the genrator

% Calculate the electric energy at every vertex point
electricEnergy = power .* (effectivePower / totalPower);

figure(3);
trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', electricEnergy);
title('RFA power distribution at every point of mesh');


%% Calculate Temperature Distribution T

% The temperature distribution is given by the heat equation
% | dT(x,y,z,t)/dt - div ( (lambda(x,y,z,t) * grad T(x,y,z,t) ) = Q_heat |
%  -> This is a second order parabolic equation
% Problem is solved with a mixed finite element method

% Define the new boundary conditions for the heat equation
 bmesh = DefineBoundaryConditions(bedges, 'temp');
 
% Define specific pde parameters for parabolic heat equation
k_Temp = @(r,z) 1; 
q_Temp = @(r,z) 0;
intyp = 1;

nu  = nu_blood;    % perfusion coefficient, prefactor to Q_perf
rho = rho_blood;   % density
c   = c_blood;     % heat capacity    

T_body = 37 + 273.5; % body temperature in Kelvin

% Set initial temperature distribution for t = 0
% (T = T_body on the entire domain)
uh0_Temp = zeros(size(power)) + T_body;
uh_Temp  = uh0_Temp;

% Calculate the right hand side of the equation with discrete point
Q_rfa   = electricEnergy;                       % heat of electrical power
Q_perf  = nu .* rho .* c .* (uh_Temp - T_body); % heat of blood perfusion
 
Q_total = Q_rfa + Q_perf;

% Assemble FEM systems of equation for temperature distribution
[Ah_heat, Mh_heat, fh_heat] ...
    = AssembCylindricHeatEquation2D(pmesh, tmesh, k_Temp, q_Temp, Q_total, intyp);

 % Add boundary conditions
[Ah_heat, fh_heat] = AddBoundaryConditionsToFEMatrix(Ah_heat, fh_heat, pmesh, tmesh, bmesh);
 
stopTheExecutionHereBreakpoint = 0;



%% Abbreviations TODO: UPDATE
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

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
t_step  = 0.50;
t_end   = 5.00;  
t_vec   = t_start:t_step:t_end;


%% Material parameters

% For the first simulation run, everything should be absolutely const.

% Alternative model for linear dependent material parameters
% for ref. see books from Stein T. and Kröger T.
sigma_stein  = [ 0.21     0.013      -1     1.143 ]'; 
rho_stein    = [ 1080   -0.00056   -0.657     0   ]';
c_stein      = [ 3455*0.001 0      -0.376     0   ]';   % Stein gives c in J/gK not in J/kg*K
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


%% Plot the mesh, for control -> deactivated by comments
%figure(1);
%subplot(1,2,1);
%trimesh(tmesh', pmesh(1,:)', pmesh(2,:)');
%title('The used triangulation of the domain');

%subplot(1,2,2);
%scatter(pmesh(1,:), pmesh(2,:));
%title('All the points in the triangulation');


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
power_setup = 20;   % power of the generator (in range 20-200 W)
U_elec = 2;         % Potential difference of the two electrodes

R_setup = 100; % TODO find good value % inner resistance of the generator
R_tis = U_elec * U_elec / totalPower; % tissue resistance

effectivePower = (4 * power_setup * R_tis * R_setup) / ...
    ((R_tis + R_setup)^2);  % effective power of the genrator

% Calculate the electric energy at every vertex point
electricEnergy = power .* (effectivePower / totalPower);

%% Plot the power distribution - deactivated by comments
%figure(3);
%trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', electricEnergy);
%title('RFA power distribution at every point of mesh');


%% Calculate the Heat for the Temperature Distribution T

% The temperature distribution is given by the heat equation
% | dT(x,y,z,t)/dt - div ( (lambda(x,y,z,t) * grad T(x,y,z,t) ) = Q_heat |
%  -> This is a second order parabolic equation
% Problem is solved with a mixed finite element method

% Define the new boundary conditions for the heat equation
 bmesh = DefineBoundaryConditions(bedges, 'temp');
 
% Define specific pde parameters for parabolic heat equation
k_Temp = @(r,z) lambda; 
q_Temp = @(r,z) 0;
intyp = 1;

nu  = nu_blood;    % perfusion coefficient, prefactor to Q_perf
rho = rho_blood;   % density
c   = c_blood;     % heat capacity    
lam = lambda;      % thermal conductivity 

T_body = 37 + 273.5; % body temperature in Kelvin

% Set initial temperature distribution for t = 0
% (T = T_body on the entire domain)
uh0_Temp = zeros(size(power)) + T_body;
uh_Temp  = uh0_Temp;

% Calculate the right hand side of the equation with discrete point
Q_rfa   = electricEnergy;                       % heat of electrical power
Q_perf  = nu .* rho .* c .* (uh_Temp - T_body); % heat of blood perfusion
 
Q_total = 0;

%% Calculate the temperatute distribution over time

for t_count=2:size(t_vec,2)

    delta_t = t_vec(t_count) - t_vec(t_count-1); 
        
    % We will reduce the time dependency to a semi-discrete problem
    %   Mh * du/dt + Ah * u = fh 
    % 
    % This is now practically an ODE we can solve over time 


    % In theory, you could solve this ODE by inverting Mh
    % TODO 


    % To solve the system of ODE aborve, we use implicit Euler
    % 
    % Mh + delta_t * Ah * u(t_1) = u(t_0) + delta_t * f(t_1)
    %

    
    Q_perf  = nu .* rho .* c .* (uh_Temp - T_body); % heat of blood perfusion
 
    Q_total = Q_total + delta_t * (Q_rfa + Q_perf); % Update Q_total

    % Test
    % Q_total = zeros(size(Q_total));
    
    % Assemble FEM systems of equation for temperature distribution
    [Kh_heat, Mh_heat, fh_heat] ...
        = AssembCylindricHeatEquation2D(pmesh, tmesh, k_Temp, q_Temp, Q_total, intyp);

    %% TRY SOLVE PDE SYSTEM OLD
%     Mh_heat = Mh_heat * c * rho;
%  
%     Ah_Temp = Mh_heat + delta_t * Ah_heat;
%     fh_Temp = uh0_Temp + delta_t * fh_heat;
% 
%     % Add boundary conditions
%     [Ah_Temp, fh_Temp] = AddBoundaryConditionsToFEMatrix(Ah_Temp, fh_Temp, pmesh, tmesh, bmesh);
% 
%     % This is just for visualization
%     uh_Temp = Ah_Temp \ fh_Temp; 
% 
%     % use new solution as old solution for next iteration
%     uh0_Temp = uh_Temp;
%     
%     figure(4);
%     trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', uh_Temp - 273.15);
%     title('Schematic Temperature Distribution in ° Celsius');

    %% TRY SOLVE PDE NEW
    
    % Reduce this to a system of ODE 
    % u'(t) = A * u(t) + g(t)
    % u(0) = u_0
    
    % u(t) = u_h
    % A    = - M_h^(-1) * Kh
    % g(t) = M_h^(-1) * f_h(t)
    % u_0  = M_h^(-1) * Anfangswerte
    
    % Calculate inverted mass matrix
    Mh_heat = Mh_heat * c * rho;    
    Mh_heat_inv = inv(Mh_heat); 
    
    % Calculate the parts for ODE above  
    Ah_heat = - Mh_heat_inv * Kh_heat;
    gh_heat = + Mh_heat_inv * fh_heat; 
    uh_heat = uh0_Temp;
    
    gh_heat_neu = gh_heat;
    
    % Solve with sigma variant from Numerik 3
    n = size(Ah_heat,1);    
    sigma = 0.5;
    tau = 0.01;

    A_h = Ah_heat;

    Ah_lhs = zeros(n) + 1  + tau * sigma * A_h;        % Ah auf linker Seite
    Ah_rhs = zeros(n) + 1  - (tau * (1-sigma) * A_h);  % Ah auf rechter Seite

    % Korrektur an den Eckstellen von Ah auf rechter Seite
    %Ah_rhs(1,1) = 1;
    %Ah_rhs(n,n) = 1;    

    helper_links = Ah_rhs*uh0_Temp;
    helper_rechts = tau * (sigma*gh_heat_neu' + (1-sigma)*gh_heat');

    helper = helper_links + helper_rechts';

    [Ah_lhs, helper] = AddBoundaryConditionsToFEMatrix(Ah_lhs, helper, pmesh, tmesh, bmesh);

    uh_heat = Ah_lhs \ helper;
    
    % End solve with variant from Numerik 3

    figure(4);
    trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', uh_heat - 273.15);
    title('Schematic Temperature Distribution in ° Celsius');

end % for 



stopTheExecutionHereBreakpoint = 0;

%% Galerkin FEM - TRYING TODO maybe this is a bit to complicated ...
%  For every time step, there will be a system of equations 
%
%    R1 S1     u0     F1
%    R2 S2  *  u1  =  F2
%
%  which yields the solution u = u0 + u1
% 
%  with R1 = Mh + delta_t * Ah
%       S1 = Mh + 0.5 * delta_t * Ah
%       R2 = 0.5 * delta_t * Ah
%       S2 = 0.5 * Mh + (delta_t/3) * Ah
%       
%       F1 = f_rhs
%       F2 = 
%
% For reference, see 
% Diskontinuierliche GALERKIN-Verfahren 
% in Raum und Zeit zur Simulation von Transportprozessen



% Calculate the system of equations for galerkin fem
R1 = Mh + delta_t*A;
S1 = Mh + 0.5 * delta_t * A;
R2 = 0.5 * delta_t * Ah;
S2 = 0.5 * Mh + (delta_t/3) * Ah;
Ph = [R1, S1;   % full matrix
      R2, S2];
  


uh_Temp = Ph \ fh_heat;
 



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
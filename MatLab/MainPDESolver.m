
%% General information about the script

% The whole matlab code is part of an AMP application project
% Project is about numeric simulation of radio frequency ablation

% THIS SCRIPT IS NOT OPTIMIZED
% It heavily prefers intention telling over optimization

% For model reference see article of my advisor:
% Tim Kroeger et al. 
% "Numerical Simulation of Radio Frequency
%  Ablation with State Dependent Material
%  Parameters in Three Space Dimensions"

clear variables;
disp('This script runs a simulation of radio frequency ablation');


%% Explanation of computational domain

% This simulation models a needle inserted in a malignant tissue
% Due to axis symmetrie, the whole simulation is reduced to a 2D-problem
% All calculations are done in cylindrical coordinates

% Cylindrical domain for the simulation
%
%  Full cross-section       halved 
%   _____    _____           ____
%  |     |  |     |         |    |
%  |     |  |     |         |    |
%  |     |  |     |         |    |
%  |      \/      |   or   /     |
%  |              |        |     |
%  |              |        |     |
%  |______________|        |_____|
%

% Transformation reference to cylindrical coordinates 
% x = r * cos(phi);
% y = r * sin(phi);
% z = z;


%% Material parameters

% For the first simulation run, everything should be absolutely const.

% Alternative model for linear dependent material parameters
% for ref. see books from Stein T. and Kröger T.
% sigma_stein  = [ 0.21     0.013      -1     1.143 ]'; 
% rho_stein    = [ 1080   -0.00056   -0.657     0   ]';
% c_stein      = [ 3455*0.001 0      -0.376     0   ]';   % Stein gives c in J/gK not in J/kg*K
% lambda_stein = [ 0.437    0.0025      0       0   ]'; 

% Constant material parameters 
sigma_phi  =     0.21;      % electric conductivity
rho_blood  =  1080;         % density
c_blood    =  3455 * 0.001; % specific heat capacity
lambda     =     0.437;     % thermal conductivty -> in J/gK respectively in J/kg*K

nu_blood  =      0.01765;   % blood perfusion coefficient


%% Grid Generation - Choose grid for calculatione here

%   ->  Uncomment the grid you want to use for calculation
%       The variable below is for additional grid refinement

% 1.) A simple DEBUG Mesh
% [pmesh, tmesh, bedges] = GetSimpleDebugMesh();
%  numAdditionalGridRefinements = 0;

% 2.) Extra coarse full 2D cross-section
% [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraCoarse\');
%  numAdditionalGridRefinements = 2;
 
% 3.) Extra fine full 2D-cross-section (Warning, calculation is slow)
% [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraFine\');
%  numAdditionalGridRefinements = 0;

% 4.) Halved cross-section, coarse domain with prerefinement for region around electrodes
  [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Triang_Halved_Needle\');
   numAdditionalGridRefinements = 0;

   
%% Optional refinement of the initial grid

% number of refinement steps can be defined above

pmeshFiner = pmesh;
tmeshFiner = tmesh;
bedgesFiner = bedges;

for i=1:numAdditionalGridRefinements

    [pmeshFiner, tmeshFiner, bedgesFiner] =  ...
        TriangularMeshRefinement2D(pmeshFiner, tmeshFiner, bedgesFiner);

end 


%% Plot the mesh, for control -> can be deactivated by comments

%% TODO this is for presentation
% figure(55);
% %subplot(2,2,1);
% trimesh(tmesh, pmesh(:,1), pmesh(:,2));
% title('Triangulation on halved cross-section coarse');
% xlabel('r axis in meter');
% ylabel('z axis in meter');
% 
% figure(56);
% %subplot(2,2,1);
% trimesh(tmeshFiner, pmeshFiner(:,1), pmeshFiner(:,2));
% title('Triangulation on halved cross-section refined');
% xlabel('r axis in meter');
% ylabel('z axis in meter');


figure(1);
subplot(2,2,1);
trimesh(tmesh, pmesh(:,1), pmesh(:,2)');
title('Triangulation without refinement');
ylabel('Space dimensions in mm')

subplot(2,2,2);
trimesh(tmeshFiner, pmeshFiner(:,1), pmeshFiner(:,2));
title("Triangulation after " + num2str(numAdditionalGridRefinements) + " refinements");

subplot(2,2,3);
scatter(pmesh(:,1), pmesh(:,2), 8, 'filled');
title('Point vertices in the initial triangulation');

subplot(2,2,4);
scatter(pmeshFiner(:,1), pmeshFiner(:,2), 8, 'filled');
title('Point vertices in the refined triangulation');

% Rebase refined grid to be the used grid
pmesh = pmeshFiner;
tmesh = tmeshFiner;
bedges = bedgesFiner;


%% Calculate electrical potential phi 

% Define boundary conditions for the elliptical phi PDE
bmesh = DefineBoundaryConditions(bedges, 'phi');

% Problematic nodes, for ref. see function GetUndefinedBoundaryPoints(...)
undefinedNodes = GetUndefinedBoundaryPoints(bmesh); 


% Define PDE coefficients, ref. see function 'AssembCylindricLaplace2D(...)'
k_EPot  = @(r,z) 1;
q_EPot  = @(r,z) 0;
f_rhs_EPot = @(r,z) 0;
intyp = 2;

% Assemble global FEM system of equations elementwise
[Ah, fh] = AssembCylindricLaplace2D(pmesh, tmesh, k_EPot, q_EPot, f_rhs_EPot, intyp);

% Add boundary conditions
[Ah, fh] = AddBoundaryConditionsToFEMatrix(Ah, fh, pmesh, bmesh);

% Solve the system of equations
phi = Ah \ fh;


% Handle problematic bundary points by approximation of neighbour values
for i=1:size(undefinedNodes,1)
    phi(undefinedNodes(i,1)) = ...
        (phi(undefinedNodes(i,2)) + phi(undefinedNodes(i,3))) / 2;
end


figure(2);
trisurf(tmesh, pmesh(:,1), pmesh(:,2), phi);
title('Numerical FEM solution for electric potential phi'); 
xlabel('r axis');
ylabel('z axis');
colormap(jet);
colorbar('AxisLocation','in'); 
caxis([min(phi), max(phi)]);
%zlim([-1.5 1.5]);


%% Plot 3D mesh reconstruction for phi

% Domain is rotation symmetric
% We can use value for every point of 3d domain

figure(20)
[pmesh3D, uh3D, colorMap3D_Phi] = Recreate3DCylinderFromSlice(pmesh, phi, 4);
scatter3(pmesh3D(:,1), pmesh3D(:,2), pmesh3D(:,3), 5, colorMap3D_Phi, 'filled');
title('Electrical potential phi on discrete points');
xlabel('x axis in meter');
ylabel('y axis in meter');
colormap(jet);
colorbar('AxisLocation','in'); 
caxis([min(uh3D), max(uh3D)]);


%% Calculate electric energy of the electric potential

[effPowerPoints, effPowerElements] = ... 
    CalculateEffectivePower(pmesh, tmesh, bedges, phi, sigma_phi);


figure(3);
trisurf(tmesh, pmesh(:,1), pmesh(:,2), effPowerPoints);
title('Effective electric energy in VAs');
xlabel('r-axis in meter');
ylabel('z-axis in meter');
colormap(jet);
colorbar('AxisLocation','in'); 
caxis([min(effPowerPoints), max(effPowerPoints)]);


figure(600)
[pmesh3D, uh3D, colorMap3D_Energy] = Recreate3DCylinderFromSlice(pmesh, effPowerPoints, 2);
scatter3(pmesh3D(:,1), pmesh3D(:,2), pmesh3D(:,3), 5, colorMap3D_Energy, 'filled');
title('Effective electric energy');
xlabel('x axis in meter');
ylabel('y axis in meter');
colormap(jet);
colorbar('AxisLocation','in'); 
caxis([min(uh3D), max(uh3D)]);



stopExecutionHereBeforeEnteringTemperatureCalculation = 0;


%% Calculate the heat used for the temperature distribution T

% The temperature distribution is given by the heat equation
% | dT(x,y,z,t)/dt - div ( (lambda(x,y,z,t) * grad T(x,y,z,t) ) = Q_heat |
%  -> This is a second order parabolic equation
% Problem is solved with a mixed finite element method

kelvinToCelsius = 273.15;

% Define the new boundary conditions for the heat equation
 bmesh = DefineBoundaryConditions(bedges, 'temp');
 undefinedNodes = GetUndefinedBoundaryPoints(bmesh);
 
% Define specific pde parameters for parabolic heat equation
k_Temp = @(r,z) lambda; 
q_Temp = @(r,z) 0;
intyp = 1;

% Material parameters, see section above
nu  = nu_blood;    % perfusion coefficient, prefactor to Q_perf
rho = rho_blood;   % density
c   = c_blood;     % heat capacity    
lam = lambda;      % thermal conductivity 

T_body = 37 + kelvinToCelsius; % body temperature in Kelvin

% Set initial temperature distribution for t = 0
% (T = T_body on the entire domain)
uh0_Temp = zeros(size(effPowerPoints)) + T_body;
uh_Temp  = uh0_Temp;

 
Q_total = 0;


%% Time discretization -> Set step size here

% Time discretization in seconds for time-dependant simulation 
% Since problem is stiff, smaller time step leads to better results

t_start = 0.00;    % Starting point -> t = 0 seconds
t_step  = 0.2;     % delta t on uniform time steps
t_end   = 20.00;  
t_vec   = (t_start:t_step:t_end)';


%% Calculate the temperatute distribution ODE over time


t_next = t_vec(1);
uh_next = uh0_Temp;

% Time looping

for t_count=2:1:length(t_vec)
    
    uh_old = uh_next;
    
    % Get time step delta_t
    t_old  = t_next;
    t_next = t_vec(t_count);      
    delta_t = t_next - t_old;
        
    Q_rfa = delta_t .* effPowerPoints;         % sum of electricEnergy
    Q_perf  = delta_t .* nu .* rho .* c .* (T_body - uh_next); % cooling of blood perfusion
    Q_total = Q_total + Q_rfa + Q_perf;
    
    [Kh_heat, Mh_heat, fh_heat] ...
          = AssembCylindricHeatEquation2D(pmesh, tmesh, k_Temp, q_Temp, Q_total, intyp);
          
    left  = Mh_heat + delta_t * Kh_heat;
    right = Mh_heat * uh_old + delta_t * fh_heat;
    
    [left, right] = AddBoundaryConditionsToFEMatrix(left, right, pmesh, bmesh);      
     
    uh_next = left \ right;
    
    % Approximate values on problematic points from neighbours
    for i=1:size(undefinedNodes,1)
        uh_next(undefinedNodes(i,1)) = ... 
            (uh_next(undefinedNodes(i,2)) + uh_next(undefinedNodes(i,3))) / 2;
    end
  
    % Update plot for the calculation results of temperature distribution
    figure(4);
    trisurf(tmesh, pmesh(:,1), pmesh(:,2), uh_next - kelvinToCelsius);
    title(['Temperature distribution in ° Celsius after ', num2str(t_next), ' seconds']);
    colormap(jet);
    colorbar(); 
    caxis([min(uh_next)-kelvinToCelsius, max(uh_next)-kelvinToCelsius])
    
    toBeReservedForBrakePoint = 0;
        
end % for 


%% Create 3D-Data from 2D slice

figure(1000)
[pmesh3D, uh3D, colorMap3D_Temp] ... 
    = Recreate3DCylinderFromSlice(pmesh, uh_next-kelvinToCelsius, 4);
scatter3(pmesh3D(:,1), pmesh3D(:,2), pmesh3D(:,3), 5, colorMap3D_Temp, 'filled');
colormap(jet);
colorbar(); 
caxis([min(uh3D), max(uh3D)])
title(['Temperature distribution in °C after ', num2str(t_next), ' seconds']);
xlabel('x-axis in meter');
ylabel('y-axis in meter');


stopSkriptHere = 0;


%% OLD REGION SOLVE ODE -> Experimental code fragments

% This is fragmented code from calculation of temperature distribution
% Will be documented for later use in project


% -> An approach using sigma for Euler/Crank-Nicolson

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


% Q_perf  = nu .* rho .* c .* (uh_Temp - T_body); % heat of blood perfusion
% Q_total = Q_total + delta_t * (Q_rfa + Q_perf); % Update Q_total

% Test
% Q_total = zeros(size(Q_total));

% Assemble FEM systems of equation for temperature distribution
% [Kh_heat, Mh_heat, fh_heat] ...
  %   = AssembCylindricHeatEquation2D(pmesh, tmesh, k_Temp, q_Temp, Q_total, intyp);

  
  
% TRY SOLVE Matrix by matrix invert matrix 
    
% Reduce this to a system of ODE 
% u'(t) = A * u(t) + g(t)
% u(0) = u_0

% u(t) = u_h
% A    = - M_h^(-1) * Kh
% g(t) = M_h^(-1) * f_h(t)
% u_0  = M_h^(-1) * Anfangswerte

% Calculate inverted mass matrix
%     Mh_heat = Mh_heat * c * rho;    
%     Mh_heat_inv = inv(Mh_heat); 
%     
%     % Calculate the parts for ODE above  
%     Ah_heat = - Mh_heat_inv * Kh_heat;
%     gh_heat = + Mh_heat_inv * fh_heat; 
%     uh_heat = uh0_Temp;
%     
%     gh_heat_neu = gh_heat;
%     
%     % Solve with sigma variant from Numerik 3
%     n = size(Ah_heat,1);    
%     sigma = 0.5;
%     tau = 0.01;
% 
%     A_h = Ah_heat;
% 
%     Ah_lhs = zeros(n) + 1  + tau * sigma * A_h;        % Ah auf linker Seite
%     Ah_rhs = zeros(n) + 1  - (tau * (1-sigma) * A_h);  % Ah auf rechter Seite
% 
%     % Korrektur an den Eckstellen von Ah auf rechter Seite
%     %Ah_rhs(1,1) = 1;
%     %Ah_rhs(n,n) = 1;    
% 
%     helper_links = Ah_rhs*uh0_Temp;
%     helper_rechts = tau * (sigma*gh_heat_neu' + (1-sigma)*gh_heat');
% 
%     helper = helper_links + helper_rechts';
% 
%     [Ah_lhs, helper] = AddBoundaryConditionsToFEMatrix(Ah_lhs, helper, pmesh, tmesh, bmesh);
% 
%     uh_heat = Ah_lhs \ helper;
%     
%     % End solve with variant from Numerik 3
% 
%     figure(4);
%     trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', uh_heat - 273.15);
%     title('Schematic Temperature Distribution in ° Celsius');


% -> Galerkin FEM - TRYING 
% 
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
% R1 = Mh + delta_t*A;
% S1 = Mh + 0.5 * delta_t * A;
% R2 = 0.5 * delta_t * Ah;
% S2 = 0.5 * Mh + (delta_t/3) * Ah;
% Ph = [R1, S1;   % full matrix
%       R2, S2];
%   
% 
% 
% uh_Temp = Ph \ fh_heat;
 

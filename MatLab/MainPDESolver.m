
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
% for ref. see books from Stein T. and Kr�ger T.
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

%   ->  Uncomment the grid you want to investigate

% 1.) Simple DEBUG Mesh
% [pmesh, tmesh, bedges] = GetSimpleDebugMesh();
%  numAdditionalGridRefinements = 0;

% 2.) Extra coarse full 2D Cross-Section
% [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraCoarse\');
%  numAdditionalGridRefinements = 2;
 
% 3.) Extra fine full 2D-cross-section
% [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Electrodes_Triang_ExtraFine\');
%  numAdditionalGridRefinements = 0;

% 4.) Halved Cross-Scetion, coarse domain with prerefinement for region around electrodes
  [pmesh, tmesh, bedges] = ReadGridFromFile('Grid\Unstruc_Triang_Halved_Needle\');
   numAdditionalGridRefinements = 1;


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

figure(1);
subplot(2,2,1);
trimesh(tmesh, pmesh(:,1), pmesh(:,2)');
title('Triangulation without refinement');
ylabel('Space dimensions in mm')

subplot(2,2,2);
trimesh(tmeshFiner, pmeshFiner(:,1), pmeshFiner(:,2));
title("Triangulation after " + num2str(numAdditionalGridRefinements) + " refinements");

subplot(2,2,3);
scatter(pmesh(:,1), pmesh(:,2));
title('Point vertices in the initial triangulation');

subplot(2,2,4);
scatter(pmeshFiner(:,1), pmeshFiner(:,2));
title('Point vertices in the refined triangulation');

% Rebase refined grid to be the used grid
pmesh = pmeshFiner;
tmesh = tmeshFiner;
bedges = bedgesFiner;


%% Calculate electrical potential phi 

% Add boundary conditions for the elliptical problem
bmesh = DefineBoundaryConditions(bedges, 'phi');
undefinedPoints = GetUndefinedBoundaryPoints(bmesh); % ref. see function

% Define specific parameters for phi PDE

k_EPot  = @(r,z) 1;
q_EPot  = @(r,z) 0;
f_rhs_EPot = @(r,z) 0;
intyp = 1;

% Assemble FEM system of equations 
[Ah, fh] = AssembCylindricLaplace2D(pmesh, tmesh, k_EPot, q_EPot, f_rhs_EPot, intyp);

% Add boundary conditions
[Ah, fh] = AddBoundaryConditionsToFEMatrix(Ah, fh, pmesh, bmesh);

% Solve the system of equations
phi = Ah \ fh;


% Handle problematic bundary points by approximation 
for i=1:length(undefinedPoints)
    phi(undefinedPoints(i,1)) = (phi(undefinedPoints(i,2)) + phi(undefinedPoints(i,3))) / 2;
end


figure(2);
trisurf(tmesh, pmesh(:,1), pmesh(:,2), phi);
title('Solution of the finite element method'); % for phi');
xlabel('x axis');
ylabel('y-axis');
%zlim([-1.5 1.5]);


%% Plot 3D mesh reconstruction -> TODO move downwards
% Domain is rotation symmetric
% We can use value for every point of 3d domain


% uh = zeros(size(pmesh,1),1);
% [pmesh3D, uh3D] = Recreate3DCylinderFromSlice(pmesh,uh, 4);
% d = pmesh3DCylinder;
% plot3(d(:,1), d(:,2), d(:,3));


%% Calculate electric power from the electric potential

[energyPoints, energyElements] = ... 
    CalculateElectricEnergy(pmesh, tmesh, bedges, phi, sigma_phi);


figure(3);
trisurf(tmesh, pmesh(:,1), pmesh(:,2), energyPoints);
title('RFA electric energy at every point of mesh');


%% Calculate the feat used for the temperature distribution T

% The temperature distribution is given by the heat equation
% | dT(x,y,z,t)/dt - div ( (lambda(x,y,z,t) * grad T(x,y,z,t) ) = Q_heat |
%  -> This is a second order parabolic equation
% Problem is solved with a mixed finite element method

% Define the new boundary conditions for the heat equation
 bmesh = DefineBoundaryConditions(bedges, 'temp');
 
 undefinedPoints = GetUndefinedBoundaryPoints(bmesh);
 
% Define specific pde parameters for parabolic heat equation
k_Temp = @(r,z) lambda; 
q_Temp = @(r,z) 0;
intyp = 1;

% Material parameters, see section above
nu  = nu_blood;    % perfusion coefficient, prefactor to Q_perf
rho = rho_blood;   % density
c   = c_blood;     % heat capacity    
lam = lambda;      % thermal conductivity 

T_body = 37 + 273.15; % body temperature in Kelvin

% Set initial temperature distribution for t = 0
% (T = T_body on the entire domain)
uh0_Temp = zeros(size(energyPoints)) + T_body;
uh_Temp  = uh0_Temp;

 
Q_total = 0;
Q_rfa   = 0;



%% Time discretization -> Set step size here

% Time discretization in seconds for time-dependant simulation 
t_start = 0.00;    % Starting point -> t = 0 seconds
t_step  = 0.2;     % tau
t_end   = 240.00;  
t_vec   = t_start:t_step:t_end;


%% Calculate the temperatute distribution ODE over time

% Time looping

t_next = t_vec(1);
uh_next = uh0_Temp;

for t_count=2:size(t_vec,2)
    
    uh_old = uh_next;
    
    % Calculate DELTA t
    t_old  = t_next;
    t_next = t_vec(t_count);      
    delta_t = t_next - t_old;
        
    % Calculate total energy
%     if (t_count == 2)
%         Q_rfa = electricEnergy; % time independent by now  
%         Q_total = Q_rfa;
%     end
    Q_rfa = delta_t .* energyPoints;         % sum of electricEnergy
    Q_perf  = delta_t .* nu .* rho .* c .* (T_body - uh_next); % cooling of blood perfusion
    Q_total = Q_total + Q_rfa + Q_perf;
    
    [Kh_heat, Mh_heat, fh_heat] ...
          = AssembCylindricHeatEquation2D(pmesh, tmesh, k_Temp, q_Temp, Q_total, intyp);
          
    left  = Mh_heat + delta_t * Kh_heat;
    right = Mh_heat * uh_old + delta_t * fh_heat;
    
    [left, right] = AddBoundaryConditionsToFEMatrix(left, right, pmesh, bmesh);      
     
    uh_next = left \ right;
    
    % Approximate values on problematic points from neighbours
    for i=1:length(undefinedPoints)
        uh_next(undefinedPoints(i,1)) = ... 
            (uh_next(undefinedPoints(i,2)) + uh_next(undefinedPoints(i,3))) / 2;
    end
  
    % Update plot for the calculation results of temperature distribution
    figure(4);
    trisurf(tmesh, pmesh(:,1), pmesh(:,2), uh_next - 273.15);
    title(['Temperature Distribution in � Celsius after ', num2str(t_next), ' seconds']);
    
    if (t_count == 2)
        breakPointAfter1Second = 0;
    end
    



    if (t_count == 2)
       merken1 = uh_next; 
    elseif (t_count == 240)
       merken2 = uh_next; 
    elseif (t_count == 480)
       merken3 = uh_next; 
    elseif (t_count == 720)
       merken4 = uh_next; 
    elseif (t_count == 960)
      merken5 = uh_next; 
    end
    
    stopTheExecutionHereBreakpoint = 0;
    
    % uh_next = CalculateSingleTimeStep(Ah, Mh, fh_rhs, uh_old, t_next, delta_t);
    
    
end % for 


figure(500);
trisurf(tmesh, pmesh(:,1), pmesh(:,2), uh_next);
title('Difference between 1 minute and one second');

figure(5);
trisurf(tmesh, pmesh(:,2), pmesh(:,1), merken2 - merken1);
title('Difference between 1 minute and one second');

figure(6);
trisurf(tmesh, pmesh(:,2), pmesh(:,1), merken3 - merken1);
title('Difference between 2 minutes and one second');

figure(7);
trisurf(tmesh, pmesh(:,2), pmesh(:,1), merken4 - merken1);
title('Difference between 3 minutes and one second');

% figure(8);
% trisurf(tmesh', pmesh(2,:)', pmesh(1,:)', merken5 - merken1);
% title('Difference between 4 minutes and one second');

figure(9);
trisurf(tmesh, pmesh(:,2), pmesh(:,1), uh_next - merken1);
title('Difference between 5 minutes and one second');


stopTheExecutionHereBreakpoint = 0;

test11 = uh_next - merken1;
test22 = uh_next - merken4;


stopTheExecutionHereBreakpoint = 0;

%% Create 3D-Data from 2D slice

% see: https://www.mathworks.com/matlabcentral/answers/506492-how-do-i-plot
%      -animation-of-temperature-data-for-a-3d-object-with-time
% for video idea 
% 
% maxt = max(t);
% mint = min(t);
% N = 100;
% cmap = jet(N);                % generate N colors for temperature
% plot3(x,y,z)
% hold on
% for i = 1:length(x)
%     itemp = round((t(i)-mint)/(maxt-mint)*(N-1))+1;
%     h = plot3(x(i),y(i),z(i),'.','color',cmap(itemp,:));
%     pause(0.1)
%     delete(h)
% end
% hold off

% TODO




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
%     title('Schematic Temperature Distribution in � Celsius');


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
 

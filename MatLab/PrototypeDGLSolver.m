% This Skript is for experimenting with the RFA DGLs

% for reference see article of my advisor:
% Tim Kroeger et al. 
% "Numerical Simulation of Radio Frequency
%  Ablation with State Dependent Material
%  Parameters in Three Space Dimensions"

disp('MatLab ..... yeah :(')

%% Space and time dimensions

% cylindric domain
%   ______ 
%  /      \
% |   ()   | 
%  \______/

% Transformation to zylindric coordinates 
x = @(r,phi,z) r * cos(phi);
y = @(r,phi,z) r * sin(phi);
z = @(r,phi,z) z;

% time dimension
t_start = 0.00;
t_step  = 0.01;
t_end   = 5.00;   
T = [t_start t_step t_end];

%% Material parameters
% For the first run, everything should be const.
% for ref. see books from Stein T. and Kröger T.

% Constant material parameters
sigma  = [ 0.21     0.013      -1     1.143 ]'; % electric conductivity
rho    = [ 1080   -0.00056   -0.657     0   ]'; % density
c      = [ 3455       0      -0.376     0   ]'; % specific heat capacity
lambda = [ 0.437    0.0025      0       0   ]'; % thermal conductivty

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






%% PDE

% Electric energy


% fixed potential on electrodes (I x OMEGA_elec)
% phi = +/- 1  
% -- First try only, this is very arbitrarily
% -- Actually the potential is induced by a function

% inner domain (I x (OMEGA \ OMEGA_elec) )
% - (NABLA) ( (sigma * (NABLA) phi) ) = 0

% robin boundary condition (I x GAMMA)
% n * NABLA phi = ( (n * (s-x)) / ( |s-x|^2 ) ) * phi
% s = barycenter of the union of all probes

% 
% Q_rf(t,x) = p((t,x)) * (p_eff(t)/p_total(t))

Q_rf   = @(t,x) 1; % TODO
Q_perf = @(t,x) 1; % TODO

Q_pc   = @(t,x) 0; % Is not considered in first simulation 
                   % Cause its not reliably (Kroeger)

pde_phi = @(x) 1;
pde_Q = @(t,x) Q_rf(t,x) + Q_perf(t,x) + Q_pc(t,x);



% Temperature 

% TODO

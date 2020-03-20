% This Skript is for experimenting with the RFA DGLs

% for reference see article of my advisor:
% Tim Kroeger et al. 
% "Numerical Simulation of Radio Frequency
%  Ablation with State Dependent Material
%  Parameters in Three Space Dimensions"

disp('MatLab ..... yeah :(')


% Material parameters
% For the first run, everything is const.

% for ref. see books from Stein T. and Kröger T.
sigma  = [ 0.21     0.013      -1     1.143 ]';
rho    = [ 1080   -0.00056   -0.657     0   ]';
c      = [ 3455       0      -0.376     0   ]';
lambda = [ 0.437    0.0025      0       0   ]';



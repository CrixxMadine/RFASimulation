function [pmesh3D, uh3D, colorMap3D] = Recreate3DCylinderFromSlice(pmesh2D, uh2D, rotType)

%% Function and parameter description

% Recreate the 2D cross-section solution on a scattered 3D cylinder 
% The cross section is rotated around the axis r = 0

% Returns 
% 
%  pmesh3D := new mesh with 3D cartesian coordinates (x,y,z)
%  uh3D    := solution vector for 3D points
%
%  colorMap3D := a RGB colormap corresponding to the solution uh3D
%                this comes in handy if you want to plot the 3D solution

% Input arguments:
% 
%  pmesh2D := 2D coordinates of discrete nodes
%  uh2D    := solution vector for the 2D nodes 
% 
%  rotType := a value signaling the number of discrete rotation steps
%             e.g. rotType = 1 returns discrete rotations angular 90°
% 
%  rotType = 1 :   4 points  ->  90 °
%  rotType = 2 :   9 points  ->  40 ° 
%  rotType = 3 :  18 points  ->  20 °
%  rotType = 4 :  36 points  ->  10 °
%  rotType = 5 :  72 points  ->   5 °
%  rotType = 6 : 120 points  ->   3 °


%% Implementaion

% 1.) -> Define the number of discrete rotation points

% Default value if no rotation type is defined
if (nargin == 2)
    rotType = 4;  % -> Discrete points for every 10° angle
end


switch rotType
    case 1
        angleStepSize = 90; 
    case 2
        angleStepSize = 40;
    case 3 
        angleStepSize = 20;
    case 4 
        angleStepSize = 10;
    case 5
        angleStepSize =  5;
    case 6
        angleStepSize =  3;
    otherwise
        error('You entered and invalid rotType argument!')
end

lastAngle = 360 - angleStepSize;
numberOfSteps = 360 / angleStepSize;


% 2.) -> Generate a color gradient for the values of solution vector

numberOfNodes = size(pmesh2D,1);

numberOfDifferentColors = 256;   % Default, you can use other values too!
colorGradientRGB = jet(numberOfDifferentColors);

maxVal = max(uh2D); 
minVal = min(uh2D);

% Set the minimum value in uh2D to zero for referencing color gradient
uh_Reference = uh2D - minVal;
currentColorMap2D(numberOfNodes,3) = 0;

referenceFactor = (numberOfDifferentColors - 1) / (maxVal - minVal);

for count=1:numberOfNodes
    colorRow = floor(referenceFactor * uh_Reference(count) + 1);
    currentColorMap2D(count,:) = colorGradientRGB(colorRow,:);
end

colorMap3D = zeros(numberOfNodes * numberOfSteps , 3);


% 3.) Rebuild 3D cartesian coordinates from 2D cylinder coordinates

pmesh3DCylinder = zeros(numberOfNodes * numberOfSteps , 3);
uh3D = zeros(numberOfNodes * 4, 1);

for angle=0:angleStepSize:lastAngle
    
    angleVec = zeros(numberOfNodes,1) + angle;
    
    range = (angle/angleStepSize) * numberOfNodes + 1;
    pmesh3DCylinder(range:range+numberOfNodes-1, :) = [pmesh2D(:,1), angleVec, pmesh2D(:,2)];   
      
    uh3D(range:range+numberOfNodes-1) = uh2D;
   
    colorMap3D(range:range+numberOfNodes-1, :) = currentColorMap2D; 
      
end

[s,t,u] = pol2cart(pmesh3DCylinder(:,2)', pmesh3DCylinder(:,1)', pmesh3DCylinder(:,3)');

pmesh3D = [s',t',u'];

end


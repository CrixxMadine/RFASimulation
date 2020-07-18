function [pmesh3D, uh3D, colorMap3D] = Recreate3DCylinderFromSlice(pmesh, uh, rotType)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% rotType = 1 :   4 points
% rotType = 2 :   9 points
% rotType = 3 :  18 points
% rotType = 4 :  36 points
% rotType = 5 :  72 points
% rotType = 6 : 120 points


% Default value if no rotation type is defined
if (nargin == 2)
    rotType = 4;
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

%% Implementaion

numberOfNodes = size(pmesh,1);

% 1.) -> Calculate a color gradient for the solution vector

numberOfDifferentColors = 256;  % use can use other value too 
colorGradientRGB = jet(numberOfDifferentColors);

maxVal = max(uh); 
minVal = min(uh);

uh_Reference = uh - minVal;
currentColorMap(numberOfNodes,3) = 0;

referenceFactor = (numberOfDifferentColors - 1) / (maxVal - minVal);

for count=1:numberOfNodes
    colorRow = floor(referenceFactor * uh_Reference(count) + 1);
    currentColorMap(count,:) = colorGradientRGB(colorRow,:);
end

colorMap3D = zeros(numberOfNodes * numberOfSteps , 3);


% 2.) Rebuild 3D cartesian coordinates from 2D cylinder coordinates

pmesh3DCylinder = zeros(numberOfNodes * numberOfSteps , 3);
uh3D = zeros(numberOfNodes * 4, 1);

for angle=0:angleStepSize:lastAngle
    
    angleVec = zeros(numberOfNodes,1) + angle;
    
    area = (angle/angleStepSize) * numberOfNodes + 1;
    pmesh3DCylinder(area:area+numberOfNodes-1, :) = [pmesh(:,1), angleVec, pmesh(:,2)];   
      
    uh3D(area:area+numberOfNodes-1) = uh;

    
    colorMap3D(area:area+numberOfNodes-1, :) = currentColorMap; 
    
    
end

[s,t,u] = pol2cart(pmesh3DCylinder(:,2)', pmesh3DCylinder(:,1)', pmesh3DCylinder(:,3)');

%figure(90);
%plot3(s(:), t(:), u(:), '.');
% scatter3(s(:), t(:), u(:), '.');

pmesh3D = [s',t',u'];
end


function [pmesh3D, uh3D, colorMap3D] = Recreate3DCylinderFromSlice(pmesh, uh, rotType)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% rotType = 1 :  4 points
% rotType = 2 :  9 points
% rotType = 3 : 18 points
% rotType = 4 : 36 points
% rotType = 5 : 72 points
% rotType = 6 : 108 points

numberOfNodes = size(pmesh,1);

%% TSET JET
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

colorMap3D = zeros(numberOfNodes * 36 , 3);

for xxx=0:1:35
    area = (xxx) * numberOfNodes + 1;
    colorMap3D(area:area+length(uh)-1, :) = currentColorMap;    
end

% END TEST



pmesh3DCylinder = zeros(numberOfNodes * 36 , 3);

uh3D = zeros(numberOfNodes * 4, 1);

for angle=0:10:350
    
    angleVec = zeros(numberOfNodes,1) + angle;
    
    area = (angle/10) * numberOfNodes + 1;
    pmesh3DCylinder(area:area+numberOfNodes-1, :) = [pmesh(:,1), angleVec, pmesh(:,2)];
    
    uh3D(area:area+numberOfNodes-1) = uh;

end

[s,t,u] = pol2cart(pmesh3DCylinder(:,2)', pmesh3DCylinder(:,1)', pmesh3DCylinder(:,3)');

%figure(90);
%plot3(s(:), t(:), u(:), '.');
% scatter3(s(:), t(:), u(:), '.');

pmesh3D = [s',t',u'];
end


function [pmesh3D,uh3D] = Recreate3DCylinderFromSlice(pmesh, uh, rotType)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% rotType = 1 :  4 points
% rotType = 2 :  9 points
% rotType = 3 : 18 points
% rotType = 4 : 36 points
% rotType = 5 : 72 points
% rotType = 6 : 108 points

numPoints = size(pmesh,2);

pmesh3DCylinder = zeros(numPoints * 36 , 3);

for angle=0:10:350
    
    angleVec = zeros(numPoints,1) + angle;
    
    area = (angle/10) * numPoints + 1;
    pmesh3DCylinder(area:area+numPoints-1, :) = [pmesh(1,:)', angleVec, pmesh(2,:)'];

end

[s,t,u] = pol2cart(pmesh3DCylinder(:,2)', pmesh3DCylinder(:,1)', pmesh3DCylinder(:,3)');
plot3(s(:), t(:), u(:), '.');
end


function [pmesh, tmesh, bedges] = ReadGridFromFile(folderPath)

% Reads a grid from given text file
% Expects a folder path containing text files named:
% pmesh.txt, tmesh.txt, bedges.txt, binfo.txt, bconversion.txt
%
% returns pmesh, tmesh, bedges 
% 
% folderPath := location of the folder containing the upper text files
%
% Will do more generic file input as soon as needed
% For testing purposes this is enough
%
% This is how the files must be formatted: 
% pmesh.txt  : x and y coordinates of the points
% tmesh.txt  : number of p1, p2, and p3 of the triangle
% bedges.txt : number of p1 and p2 of the boundary edge
% binfo.txt  : number of boundary or the boundary in bedges
% bconversion.txt : describes which boundary number represents what domain
%    domain1   -> '100' : positive electrode
%    domain2   -> '101' : negative electrode
%    domain3   -> '200' : blank needle
%    domain4   -> '300' : outer boundary
%    domain5   -> '400' : rotation axis


pmesh = load(strcat(folderPath, 'pmesh.txt'));

tmesh = load(strcat(folderPath, 'tmesh.txt'));

% MatLab counts indizes from one
% Grid Data counts indizes from zero
tmesh  = ShiftAllPointsByOne(tmesh);


bpoints = load(strcat(folderPath, 'bedges.txt'));
bpoints = ShiftAllPointsByOne(bpoints);

bnumbers = load(strcat(folderPath, 'binfo.txt'));
bnumbers = ShiftAllPointsByOne(bnumbers);

bconverter = load(strcat(folderPath, 'bconversion.txt'));

if (max(bnumbers) > size(bconverter))
    throw(MException(folderPath,'Conversion data is corrupted'));
end

btype = ConvertBoundaryType(bnumbers, bconverter(:,2))';

bedges = [bpoints btype]; 

%% Helper functions

function matrix = ShiftAllPointsByOne(matrix)        
    matrix = matrix + 1;
end

function newInfo = ConvertBoundaryType(boundaries, conversion)
    n = size(boundaries);
    newInfo(n) = 0;

    for i=1:n
        number = boundaries(i);
        newInfo(i) = conversion(number); 
    end

end

end


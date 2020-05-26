function [pmesh, tmesh, bedges] = ReadGridFromFile(folderPath)

% Reads a grid from given text file
% Expects a folder path containing text files named:
% pmesh.txt, tmesh.txt, bedges.txt, binfo.txt
%
% returns pmesh, tmesh, bedges 
% 
% folderPath := location of the folder containing the upper text files
%
% Will do more generic file input as soon as needed
% For testing purposes this is enough

tmesh = load(append(folderPath, 'tmesh.txt'));
tmesh = int32(tmesh);

pmesh = load(append(folderPath, 'pmesh.txt'));

temp1 = load(append(folderPath, 'bedges.txt'));
temp2 = load(append(folderPath, 'binfo.txt'));

bedges = int32([temp1 temp2]);

end


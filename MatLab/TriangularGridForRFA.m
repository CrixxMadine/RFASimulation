function [pmesh, tmesh, edges] = TriangularGridForRFA()

%% Info on this subroutine
%
% This class was used to create a structured triangular grid
% Meanwhile, simulation uses unstructured grid
% Function is obsolete now

%% Parameters Explanation
%  pmesh := Point mesh, the x and y coords to all nodes 
%           e.g. pmesh(1) = [x y]
%
%  tmesh := Triangle mesh, the 3 global nodes of every triangle
%           nodes counted counterclockwise
%           e.g. tmesh(1) = [P1 P2 P3]
%           
%  edges := edges and relvant info on boundary conditions
%           e.g. edges(1) = [P1 P2 b]
%           b = 0 : inner edge
%           b = 1 : edge that touches needle
%           b = 2 : edge that touches outer boundary

%%  Structured grid variant explanation
%
%   Grid Variant 1, triangular structured
%   _______________     _______________
%  | 1/|3 /|5 /|7 /|   | 9/|11/|13/|15/|     
%  | /2| /4| /6| /8|   | / | / | / | / |   
%  |/__|/__|/__|/__|   |/__|/__|/__|/__|  
%  |17/|19/|21/|23/|   |25/|27/|29/|31/|     
%  | / | / | / | / |   | / | / | / | / |   
%  |/__|/__|/__|/__|   |/__|/__|/__|/__|  
%  |  /|  /|  /|  /|   |  /|  /|  /|  /|
%  | / | / | / | / |   | / | / | / | / |      
%  |/__|/__|/__|/__|___|/__|/__|/__|/__|
%  |49/|  /|53/|  /|57/|  /|61/|  /|65/|
%  | / | / | / | / | / | / | / | / | / |
%  |/__|/__|/__|/__|/__|/__|/__|/__|/__|       
%  |67/|  /|71/|  /|75/|  /|79/|  /|83/|
%  | / | / | / | / | / | / | / | / | / |
%  |/__|/__|/__|/__|/__|/__|/__|/__|/__| 
%  |85/|  /|89/|  /|93/|  /|97/|99/|  /|
%  | / | / | / | / | / | / | / | / | / |
%  |/__|/__|/__|/__|/__|/__|/__|/__|/__|  


%   Grid Variant 2, squared and structured
%   ______ ______ ______ ______        ______ ______ ______ ______
%  |      |      |      |      |      |      |      |      |      |
%  |   1  |   2  |   3  |   4  |      |   5  |   6  |   7  |   8  |
%  |______|______|______|______|   N  |______|______|______|______|
%  |      |      |      |      |      |      |      |      |      |
%  |   9  |  10  |  11  |  12  |   E  |  13  |  14  |  15  |  16  |
%  |______|______|______|______|      |______|______|______|______|
%  |      |      |      |      |   E  |      |      |      |      |
%  |  17  |      |  19  |      |      |  21  |      |  23  |      |
%  |______|______|______|______|   D  |______|______|______|______|
%  |      |      |      |      |      |      |      |      |      |
%  |  25  |      |  27  |      |   L  |  29  |      |  31  |      |
%  |______|______|______|______|      |______|______|______|______|
%  |      |      |      |      |   E  |      |      |      |      |
%  |  33  |      |  35  |      |      |  37  |      |  39  |      |
%  |______|______|______|______|______|______|______|______|______|
%  |      |      |      |      |      |      |      |      |      |
%  |  41  |      |  43  |      |  45  |      |  47  |      |  49  |
%  |______|______|______|______|______|______|______|______|______|
%  |      |      |      |      |      |      |      |      |      |
%  |  50  |      |  52  |      |  54  |      |  56  |      |  58  |
%  |______|______|______|______|______|______|______|______|______|
%  |      |      |      |      |      |      |      |      |      |
%  |  59  |      |  61  |      |  63  |      |  65  |      |  67  |
%  |______|______|______|______|______|______|______|______|______|
%  |      |      |      |      |      |      |      |      |      |
%  |  68  |      |  70  |      |  72  |      |  74  |      |  76  |
%  |______|______|______|______|______|______|______|______|______|


%   Grid Variant 3, split the squares into triangles as shown below
%   Generate uniform triangle grid from the quadratic grid
%   Use Grid Variant 2 as starting Grid
%   _______________
%  |\             /|
%  |   \   2   /   |
%  |      \ /      |
%  |  3    x    1  |
%  |      / \      |
%  |   /   4   \   |
%  |/_____________\|


%% Grid generation

% This section was moved elsewhere

pmesh = 0;
tmesh = 0;
edges = 0;

end


function bmesh = DefineBoundaryConditions(bedges, type)

% Define the boundary conditions for the current mesh grid
% 
% returns
% bmesh -> third row is type (1=Dirichlet, 2=Neumann)
%       -> fourth row is corresponding value
% 
% input args:
% type: the type of the boundary equation to be added
% -> 'phi'  : add conditions for electrode 
% -> 'temp' : add boundary conditions for heat equation
% -> more will be added within time!

% bedges: 
% Third row in bedges explained:
% Every boundary domain is represented by a number
%    domain1   -> '100' : positive electrode
%    domain2   -> '101' : negative electrode
%    domain3   -> '200' : blank needle
%    domain4   -> '300' : outer boundary
% 
%  TODO : will define robin and needle in time

n = size(bedges, 2);

bmesh(4,n) = 0;
bmesh(1,:) = bedges(1,:);
bmesh(2,:) = bedges(2,:);

for i=1:n
    
    if (strcmp(type, 'phi')) % Add boundary vlaues for phi 
    
        if (bedges(3,i) == 100) % positive electrode
            bmesh(3,i) = 1;     % Dirichlet
            bmesh(4,i) = 1;     % potential in Voltage

        elseif (bedges(3,i) == 101) % negative electrode
            bmesh(3,i) =  1;
            bmesh(4,i) = -1;

        else % blank needle or outer boundary
            bmesh(3,i) = 2;     % Neumann 
            bmesh(4,i) = 0;
            
            % Testing
%             bmesh(3,i) = 1;     % Dirich
%             bmesh(4,i) = 0;
        end
    
    elseif strcmp(type, 'temp')
        
        %if (2 < 1)  %% TODO TESTING
        if (bedges(3,i) == 300) % outer boundary
            bmesh(3,i) = 2;     % Neumann 
            bmesh(4,i) = 0;
            
        else     % anywhere on the needle, blank or electrode
            bmesh(3,i) = 1;           % Dirichlet
            bmesh(4,i) = 37 + 273.15; % body temperature in Kelvin
        end
        
    else   
        
        error('You did not specify an allowed equation type!')
        
    end % if type
     
end % for
end


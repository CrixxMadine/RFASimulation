function bmesh = DefineBoundaryConditions(bedges, type)

%% Function and parameter description

% Define the boundary conditions for the current mesh grid
% Adds third and fourth row with info on bmesh
 
% returns
%  bmesh -> third row is type (1=Dirichlet, 2=Neumann)
%        -> fourth row is corresponding value
 
% input args:
%  type: the type of the boundary equation to be added
%  -> 'phi'  : add conditions for electrode 
%  -> 'temp' : add boundary conditions for heat equation

% bedges: 
%  Third row in bedges explained:
%  Every boundary domain is represented by a number
%    domain1   -> '100' : positive electrode
%    domain2   -> '101' : negative electrode
%    domain3   -> '200' : blank needle
%    domain4   -> '300' : outer boundary
%    domain5   -> '400' : rotation axis


%% Implementation

n = size(bedges, 1);

bmesh(n,4) = 0;
bmesh(:,1) = bedges(:,1);
bmesh(:,2) = bedges(:,2);

for i=1:n
    
    if (strcmp(type, 'phi')) % Add boundary vlaues for phi 
    
        if (bedges(i,3) == 100) % positive electrode
            bmesh(i,3) = 1;     % Dirichlet
            bmesh(i,4) = 1;     % potential in Voltage

        elseif (bedges(i,3) == 101) % negative electrode
            bmesh(i,3) =  1;
            bmesh(i,4) = -1;

        % Optional -> Simulate mass to nullify potential   
       % elseif (bedges(i,3) == 300) % Outer boundary              
       %      bmesh(3,i) = 1;        % Dirich
       %      bmesh(4,i) = 0;
        
        else % blank needle or outer boundary or rotation axis
            bmesh(i,3) = 2;     % Neumann 
            bmesh(i,4) = 0;
            

        end
    
    elseif strcmp(type, 'temp')
        
        if (bedges(i,3) == 300) % outer boundary
            bmesh(i,3) = 2;     % Neumann 
            bmesh(i,4) = 0;
         
        elseif (bedges(i,3) == 400) % rotation axis
            bmesh(i,3) = 2;         % Neumann 
            bmesh(i,4) = 0;
         
       % Variant one, probe is not cooled    
       % else                   % anywhere else
       %     bmesh(i,3) = 2;    % Neumann 
       %     bmesh(i,4) = 0;            
       
       % Variant two, cooled probe
        else     % anywhere on the needle, blank or electrode
            bmesh(i,3) = 1;           % Dirichlet
            bmesh(i,4) = 37 + 273.15; % body temperature in Kelvin
        end
        
    else   
        
        error('You did not specify an allowed type! Use phi or temp')
        
    end 
     
end 

end


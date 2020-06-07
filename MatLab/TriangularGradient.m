function [df_x,df_y] = TriangularGradient(tmeshInput , pmesh , z_val , optSurface)

%% Function summary and arguments description
%
%  This function calculates the surface gradient or the vertex gradient
%  for a given triangulation. Surface gradients are necessary to compute 
%  the vertex gradients, so this function can calculate both.
%  The gradient of each vertex is calculated by weighting the gradients
%  of the neighboring triangles in reach of the inverse distance
%  from the vertex to the center of each triangle
%
%  returns
%  dfx, dfy := numerical values of surface or vertex gradient  
%
%  input paramaters
%  tmesh := triangle matrix of triangulation
%  pmesh := point matrix of triangulation
%  z_val := discrete values for every point of the triangulation

%  optSurface := optional argument to only calculate the surface gradient
%                write 'surface' if you only want the surface gradient
%  gradients
%  For reference: The whole algorithm is out from a book
%    Duane Hanselman, Bruce Uttlefield
%    Department of Electrical and Computer Engineering
%    University of Maine 
%    'Mastering MatLab 7'
%    Pearson, 01.11.2004, First Edition


% Check if optional argument was used
if nargin==3      % you must only specify if you want surface gradient
   optSurface='vertex';
elseif nargin<3
	error('You did not provide enough input arguments.')
end
if ~ischar(optSurface)
   error('Last input argument must be a string.')
end


%% Calculate the surface gradient 

% simplify argument names
tmesh = tmeshInput';
x = pmesh(1,:)'; 
y = pmesh(2,:)';
z = z_val;

vertexNumber = length(x);

% Check for valid input data
if ~isequal(vertexNumber, length(z))
   error('Your points does not match the number of provided z values')
end

if size(tmesh,2)~=3 || any(tmesh(:)<0) || any(tmesh(:)> vertexNumber)
   error('Your triangulation does not match the point mesh length')
end

t  = tmesh(:,[1 2 3 1]);
dy = diff(y(t),1,2);                             % [y2-y1 y3-y2 y1-y3]
dx = diff(x(t),1,2);                             % [x2-x1 x3-x2 x1-x3]
delta = -sum(dy(:,[2 3 1]).*x(tmesh),2);         % calculate determinants
dx_surface = -sum(dy(:,[2 3 1]) .* z(tmesh),2) ./ delta;  % surface gradient x
dy_surface =  sum(dx(:,[2 3 1]) .* z(tmesh),2) ./ delta;  % surface gradient y

if strncmpi(optSurface,'surface',min(4,length(optSurface)))
   df_x = dx_surface;
   df_y = dy_surface;
   
   return    % vertex gradient calculation is not required
end


%% Calculate vertex gradients

df_x = zeros(vertexNumber,1);   % allocate space for results
df_y = zeros(vertexNumber,1);

xCenter = sum(x(tmesh),2)/3;   % center points of the triangles
yCenter = sum(y(tmesh),2)/3;

[vert,idx]=sort(tmesh(:));     % sort vertices in ascending order

triCount = size(tmesh,1);           
cornerNumber = [1:triCount 1:triCount 1:triCount]'; % triangle number for each vertex
cornerNumber = cornerNumber(idx);                   % shuffle triangle number to match vertices
lastVertex   = find([diff(vert);1]);                % index of last vertex element in triCount
firstVertex  = 1;

% Calculate gradient for every vertex

for counter=1:vertexNumber	                     
    
   % triangles having the vertex at the current counter
   vertexHelper = cornerNumber(firstVertex:lastVertex(counter));
   
   % starting index for next vertex
   firstVertex = lastVertex(counter) + 1;                      
   
   % inverse distances from vertex
   idist = 1 ./ hypot(x(counter) - xCenter(vertexHelper), y(counter) - yCenter(vertexHelper)); 
   
   % inverse distance weightings
   weightings=idist.'./sum(idist);        
   
   df_x(counter) = weightings * dx_surface(vertexHelper);
   df_y(counter) = weightings * dy_surface(vertexHelper);
end

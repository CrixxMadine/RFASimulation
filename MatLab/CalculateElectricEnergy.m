function electricEnergy = CalculateElectricEnergy(pmesh, tmesh, bedges, phi, sigma)

%% Calculate electric power from the electric potential

power = zeros(size(phi,1),1);

% Get the numerical gradient of every vertex
%[phi_dx, phi_dy] = TriangularGradient(tmesh, pmesh, phi);

[phi_dx, phi_dy] = TriGradient(pmesh(:,1),pmesh(:,2), phi);

% Eliminate gradient on boundary nodes -> is not defined in weak form
boundaryNodes = unique([bedges(1,:), bedges(2,:)]);
phi_dx(boundaryNodes) = 0;
phi_dy(boundaryNodes) = 0;

% Calclulate power(r,z) for every vertex

for i=1:size(power,1)
    power(i) = sigma * norm([phi_dx(i), phi_dy(i)])^2;   
end

% TODO testing
figure(200);
trisurf(tmesh, pmesh(:,1), pmesh(:,2), power);
title('Constant power at every point of mesh');

% Calculate total power of the domain
totalPower = SurfaceIntegralTriangles(tmesh, pmesh, power);

% Calculate effective power of the model 
power_setup = 200;   % power of the generator (in range 20-200 W)
U_elec = 2;          % Potential difference of the two electrodes

R_setup = 80;        % TODO find good value % inner resistance of the generator
R_tis = U_elec * U_elec / totalPower; % tissue resistance

effectivePower = (4 * power_setup * R_tis * R_setup) / ...
    ((R_tis + R_setup)^2);  % effective power of the genrator

% Calculate the electric energy at every vertex point
electricEnergy = power .* (effectivePower / totalPower);

end


function electricEnergy = CalculateElectricEnergy(pmesh, tmesh, phi, sigma)

%% Calculate electric power from the electric potential

power = zeros(size(phi,1),1);

% Get the numerical gradient of every vertex
[phi_dx, phi_dy] = TriangularGradient(tmesh, pmesh, phi);

% Calclulate power(r,z) for every vertex

for i=1:size(power,1)
    power(i) = sigma * norm([phi_dx(i), phi_dy(i)])^2;   
end

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


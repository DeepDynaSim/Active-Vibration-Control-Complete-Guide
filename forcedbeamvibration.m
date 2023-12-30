% Adjusted Parameters 
E = 30e9; % Young's modulus (Pa) - Reduced
I = 1e-8; % Second moment of inertia (m^4) - Reduced
rho = 7800; % Density (kg/m^3)
S = 0.01; % Cross-sectional area (m^2)
L = 2; % Length of the beam (m)
gamma = 0.01; % Damping coefficient (1/s)
z0 = 0.01; % Amplitude of oscillation at z=0 (m)
w = 2 * pi * 10; % Frequency of oscillation (rad/s)

% Discretization
Nz = 100; % Number of spatial steps
Nt = 1000; % Number of time steps
dz = L / Nz;
dt = 0.01;
z = linspace(0, L, Nz);
t = linspace(0, Nt*dt, Nt);

% Initial conditions
Fz = sin(pi * z / L); % Example initial deflection
Gz = zeros(size(z)); % Initial velocity

% Initialize u
u = zeros(Nz, Nt);

% Apply initial conditions
u(:, 1) = Fz;
u(2:end-1, 2) = Fz(2:end-1) + dt * Gz(2:end-1); % Forward Euler step for first time step

% Time-stepping
for n = 2:Nt-1
    for j = 3:Nz-2 % Adjusted to avoid index out of bounds
        % Central difference in space, Forward Euler in time
        uxx = (u(j+1, n) - 2*u(j, n) + u(j-1, n)) / dz^2;
        uxxxx = (u(j+2, n) - 4*u(j+1, n) + 6*u(j, n) - 4*u(j-1, n) + u(j-2, n)) / dz^4;
        utt = (u(j, n+1) - 2*u(j, n) + u(j, n-1)) / dt^2;
        u(j, n+1) = (rho * S / E * I) * dt^2 * (uxxxx - gamma * utt) + 2*u(j, n) - u(j, n-1);
    end
    
    % Boundary conditions
    u(1, n+1) = z0 * exp(1i * w * t(n+1)); % Driven end
    u(end, n+1) = 0; % Free end
end

% Plotting the results
mesh(t, z, real(u)); % Use real part of u
xlabel('Time (s)');
ylabel('Position along beam (m)');
zlabel('Displacement (m)');
title('Vibration of a Cantilever Beam');

% Calculating c
c = E*I / (rho*S);

% Function to find roots of the transcendental equation
Delta = @(q) 2 * (1 + cosh(q*L).*cos(q*L));

% Range for q to find roots (this range may need adjustments based on specific problem parameters)
qRange = linspace(0, 100, 10000);

% Finding roots of the transcendental equation
qRoots = [];
for i = 1:length(qRange)-1
    if sign(Delta(qRange(i))) ~= sign(Delta(qRange(i+1)))
        root = fzero(Delta, [qRange(i), qRange(i+1)]);
        qRoots = [qRoots, root];
    end
end

% Calculate eigenfrequencies
eigenfrequencies = [];
for i = 1:length(qRoots)
    q = qRoots(i);
    w = sqrt((q^4 * c^2) / (1 - 1i * gamma / w));
    eigenfrequencies = [eigenfrequencies, w];
end

% Display the eigenfrequencies
disp('Eigenfrequencies (rad/s):');
disp(eigenfrequencies);
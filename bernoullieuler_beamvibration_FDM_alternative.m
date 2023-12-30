% Parameters
L = 41.4e-2;         % Beam length in meters
wd = 14.88e-3;       % Beam width in meters
t = 1e-3;            % Beam thickness in meters
E = 70e9;            % Modulus of elasticity in N/m^2
ro = 5186.5;         % Linear mass density in kg/m
I = (wd * t^3) / 12; % Moment of inertia in m^4
A = wd * t;          % Cross-sectional area in m^2
nx = 15;             % Length of space domain (from example)
nt = 10^5;           % Length of time domain (from example)
tf = 15;             % Time of simulation (from example)

% Compute the mesh spacing and time step
dx = L / nx; 
dt = tf / nt;

% Drawing coordinates
x = linspace(0, L, nx + 1); 
t = linspace(0, tf, nt + 1);

% Create memory to save data
w = zeros(nx + 1, nt + 1); 

% Initial Conditions from example code
P = 5 * 9.81; % Tip payload
for i = 1:nx + 1
    w(i, 1) = -(P * L / (E * I)) * ((x(i)^2) / 2 - (x(i)^3) / (6 * L));
    w(i, 2) = w(i, 1);
end

% Damping, payload and other parameters (from example)
ksi = 0.2; % Damping ratio
Ms = 5; % Mass of the tip payload
alpha = (dt * 2 * ksi * sqrt(E * I / (A * ro)) / (ro * A));

% Boundary disturbance (from example)
d = 0.1 + 0.1 * sin(pi * (0:nt) * dt) + 0.1 * sin(2 * pi * (0:nt) * dt) + 0.1 * sin(3 * pi * (0:nt) * dt);

% Main loop for FDM (adapted from example)
for j = 3:nt + 1     
    for i = 3:nx - 1
        wxxxx = (w(i + 2, j - 1) - 4 * w(i + 1, j - 1) + 6 * w(i, j - 1) - 4 * w(i - 1, j - 1) + w(i - 2, j - 1)) / dx^4;
        w(i, j) = ((2 + alpha) * w(i, j - 1) - w(i, j - 2) + (-E * I * wxxxx) * (dt^2 / (ro * A))) / (1 + alpha);
    end

    % Boundary Conditions (from example)
    wxL = (w(nx + 1, j - 1) - w(nx, j - 1)) / dx;
    wxxxL = (w(nx + 1, j - 1) - 3 * w(nx, j - 1) + 3 * w(nx - 1, j - 1) - w(nx - 2, j - 1)) / dx^3;
    w(nx + 1, j) = 2 * w(nx + 1, j - 1) - w(nx + 1, j - 2) + (d(j - 1) + E * I * wxxxL) * dt^2 / Ms;
    w(nx, j) = (w(nx + 1, j) + w(nx - 1, j)) / 2; 
end

% Visualization (Mesh plot)
figure(1)
mesh(x, t, w');
title('Displacement of Beam with Boundary and Initial Conditions');
ylabel('Time (s)', 'FontSize', 12);
xlabel('x (m)', 'FontSize', 12);
zlabel('Displacement (m)', 'FontSize', 12);
view([60 45]);

% Visualization (Animation)
figure;
for j = 1:length(t)
    plot(x, w(:, j));
    xlabel('x (m)');
    ylabel('Displacement (m)');
    title(sprintf('Beam Displacement at Time = %.2f s', t(j)));
    axis([0 L min(w(:)) max(w(:))]);
    drawnow;
    pause(0.01); % Pause for a short time to slow down the animation
end
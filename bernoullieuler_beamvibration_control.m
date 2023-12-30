% Bernoulli-Euler Beam Vibration with FDM and Adjusted Robust Boundary Control for Stability

% Parameters
L = 41.4*10^-2;  % Beam length in meters
wd = 14.88*10^-3; % Beam width in meters
t = 1*10^-3;     % Beam thickness in meters
E = 70*10^9;     % Modulus of elasticity in N/m^2
ro = 5186.5;     % Linear mass density in kg/m^2
I = ((wd*(t^3))/12); % Moment of inertia in kgm^2
A = wd*t;        % Cross-sectional area in m^2
lamda = A*ro; 
ksi = 0.2;       % Damping ratio
P = 0.05;      % Tip payload in N
nt = 2*10^5;       % Length of time domain
nx = 15;         % Length of space domain
tf = 1;         % Time of simulation in seconds
g = 9.81;

% Adjusted Control Gains for Stability
k1 = 0.01;  % Reduced proportional control gain
k2 = 0.02;  % Reduced derivative control gain
k = 100;    % Reduced feedback gain

% Compute the mesh spacing and time step
dx = L/nx;       % Space step
dt = tf/nt;      % Time step

% Create mesh for drawing
x = linspace(0, L, nx+1); % Space mesh
time = linspace(0, tf, nt+1); % Time mesh

% Preallocate matrices for efficiency
w = zeros(nx+1, nt+1); % Displacement
d = zeros(nt+1, 1);   % Boundary disturbance
f = zeros(nx+1, nt+1); % External disturbance

% Natural frequencies and damping coefficient
k_vals = [1.875104069, 4.694091133, 7.854757438, 10.99554073, 14.13716839, 17.27875953] / L;
w_n = sqrt((k_vals.^4) * E * I / lamda);
ra = (2 * ksi * w_n(1) * ro * A);

% Boundary and External disturbance definitions
for j = 1:nt+1
    d(j) = 0.1 + 0.1*sin(pi*(j-1)*dt) + 0.1*sin(2*pi*(j-1)*dt) + 0.1*sin(3*pi*(j-1)*dt);
end

for i = 1:nx+1
    for j = 1:nt+1
        f(i,j) = (1 + sin(0.1*pi*(i-1)*dx*(j-1)*dt) + sin(0.2*pi*(i-1)*dx*(j-1)*dt) + sin(0.3*pi*(i-1)*dx*(j-1)*dt)) * (i-1)*dx / (20*L);
    end
end

% Initial Conditions
scaling_factor = 0.5; % Adjust this factor to scale down the initial displacement
for i = 1:nx+1
    initial_deflection = -(P*L/(E*I)) * (((x(i).^2)./2) - ((x(i).^3)./(6*L)));
    w(i,1) = scaling_factor * initial_deflection;
    w(i,2) = w(i,1); % Assuming zero initial velocity
end

% Damping coefficient for time integration
%alpha = (dt * ra / (ro * A));
alpha=0;

% Constants for boundary conditions
T = 2;  % Tension

% Control point
%control_point = round(nx / 2);  % Middle of the beam
control_point=3;

% Main loop for finite difference method

for j = 3:nt+1     
    for i = 3:nx-1
        wxxxx = (w(i+2,j-1) - 4*w(i+1,j-1) + 6*w(i,j-1) - 4*w(i-1,j-1) + w(i-2,j-1)) / dx^4;
        w(i,j) = ((2 + alpha) * w(i,j-1) - w(i,j-2) + (-E*I*wxxxx + f(i,j)) * (dt^2 / (ro * A))) / (1 + alpha);
    end
    w(nx+1,j) = 2*w(nx,j) - w(nx-1,j);
    
    % Apply control law at the control point

    wx_control = (w(control_point+1,j-1) - w(control_point,j-1)) / dx;
    wxxx_control = (w(control_point+1,j-1) - 3*w(control_point,j-1) + 3*w(control_point-1,j-1) - w(control_point-2,j-1)) / dx^3;
    wt_control = (w(control_point+1,j-1) - w(control_point+1,j-2)) / dt;
    wtx_control = (wt_control - (w(control_point,j-1) - w(control_point,j-2)) / dt) / dx;
    wtxxx_control = (wtx_control - ((w(control_point-1,j-1) - w(control_point-1,j-2)) / dt - (w(control_point-2,j-1) - w(control_point-2,j-2)) / dt) / dx^2) / dx;
    u = -E*I*wxxx_control + T*wx_control - (P/g) * (k1*wtx_control - k2*wtxxx_control) - k*(wt_control - k2*wxxx_control + k1*wx_control) - d(j-1);
    w(control_point,j) = w(control_point,j-1) + dt * u; % Control action applied at the control point

    % Fixed end boundary conditions at the starting point
    w(1,j) = 0; % Displacement at the fixed end is zero
    w(2,j) = 0; % Slope at the fixed end is zero   
end

% Plotting the displacement
figure(1)
mesh(x, time, w');
title('Displacement of Beam with Adjusted Control');
ylabel('Time (s)', 'FontSize', 12);
xlabel('x (m)', 'FontSize', 12);
zlabel('w(x,t) (m)', 'FontSize', 12);
view([60 45]);

% Set up the figure for animation
h = figure;
axis tight manual; % Ensure that getframe() returns a consistent size
filename = 'beam_vibration_control_animation.avi'; % Name of the video file
v = VideoWriter(filename); % Create a VideoWriter object
open(v);

frame_skip = 500; % Capture every 500th frame for the animation
short_pause = 0.1; % Short pause time, adjust as needed

% Animation of displacement with frame skipping and reduced pause time
for j = 1:frame_skip:length(time) % Loop with frame skipping
    plot(x, w(:,j));
    xlabel('x (m)');
    ylabel('Displacement (m)');
    title(sprintf('Beam Displacement at Time = %.2f s', time(j)));
    axis([0 L min(w(:)) max(w(:))]);
    drawnow;
    pause(short_pause); % Short pause for animation speed control
    % Capture the plot as a frame and write to the video
    frame = getframe(h); % Capture the current plot
    writeVideo(v, frame);
end
% Close the video file
close(v);
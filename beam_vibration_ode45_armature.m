function beam_vibration_ode45_armature
    % Parameters for Cantilever Beam Vibration Problem
    E = 30e9;      % Elastic Modulus (N/m^2)
    rho = 2700;    % Mass density (kg/m^3)
    L = 2;       % Length (m)
    b = 0.03;      % Width of the beam (m)
    h = 0.005;     % Depth of the beam (m)
    I = (b * h^3) / 12;  % Area moment of inertia (m^4)
    A = b * h;     % Cross-sectional area (m^2)

    % Mass, stiffness, and damping matrices with fixed-end boundary conditions
    M = rho * A * L / 420 * [156, 22*L, 54, -13*L; 22*L, 4*L^2, 13*L, -3*L^2; 54, 13*L, 156, -22*L; -13*L, -3*L^2, -22*L, 4*L^2];
    K = [12*E*I/L^3, 0, 0, 0; 0, 4*E*I/L, 0, 0; 0, 0, 12*E*I/L^3, 0; 0, 0, 0, 4*E*I/L]; % Adjusted stiffness matrix
    C = 0.01 * (M + K);  % Rayleigh damping

    % Time span
    t_span = [0 10];

    % Initial conditions (zero displacement and velocity)
    initial_conditions = zeros(8, 1);

    % Additional options for ode45 to reduce time step
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

    % Solve ODEs using ode45 with reduced time step
    [T, Y] = ode45(@(t, y) beam_ode_impulse(t, y, M, C, K, L), t_span, initial_conditions, options);

    % Enforce fixed-end boundary condition by setting displacement at the fixed end to zero
    for i = 1:length(T)
        Y(i, 1) = 0; % Set displacement at the fixed end to zero
    end

    % Plot Deflection vs. Time
    figure(1);
    plot(T, Y(:, 4)); % Assuming you're interested in the tip deflection
    xlabel('Time (s)');
    ylabel('Deflection (m)');
    title('Deflection vs. Time at Beam Tip');

    % Animation of the beam response
    v = VideoWriter('beam_animation.avi');
    open(v);
    figure(2);
    for i = 1:length(T)
        plot_beam(Y(i, 1:4), L);  % Function to plot the beam
        drawnow;
        frame = getframe(gcf);  % Capture the plot as a frame
        writeVideo(v, frame);  % Write the frame to the video
    end
    close(v);
end

function dydt = beam_ode_impulse(t, y, M, C, K, L)
    % State variables
    displacement = y(1:4);
    velocity = y(5:8);

    % Impulse load at the tip
    impulse_magnitude = 25;  % High magnitude for impulse
    impulse_duration = 0.01;   % Very short duration
    F_tip = impulse_magnitude * (t < impulse_duration);  % Impulse load
    F = [0; 0; F_tip; 0]; 

    % Second-order ODEs converted to first-order system
    dydt = [velocity; M \ (F - C * velocity - K * displacement)];
end

function plot_beam(displacement, L)
    scale_factor = 10; % Scale up the displacement for visualization
    x = linspace(0, L, 100);
    y = (displacement(1) + displacement(2) * x) * scale_factor;  % Apply scaling
    plot(x, y, 'b-', 'LineWidth', 2);
    xlim([0 L]);
    ylim([-5, 5]); % Adjusted for scaled displacement
    xlabel('x (m)');
    ylabel('Scaled w(x, t) (m)');
    title('Scaled Response of the Cantilever Beam to Impulse Load');
    grid on;
end
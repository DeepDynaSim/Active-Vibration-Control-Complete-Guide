function allen_cahn_1D_spectral_crank_nicolson_optimized
    % Grid and Initial Data
    N = 512; % Increase the number of spatial points for accuracy
    L = 2 * pi;
    x = linspace(0, L, N); % x as a row vector
    v = 0.25 * sin(x)'; % Define v as a column vector for proper dimensions

    epsilon = 0.001;

    % Time interval for integration
    tspan = [0, 5];
    num_frames = 200; % Increase the number of frames for smoother animation
    dt = (tspan(2) - tspan(1)) / num_frames; % Time step

    % Initialize damping coefficients (e.g., exponential decay)
    damping_coefficient_left = exp(-0.1 * x); % Damping on the left side
    damping_coefficient_right = exp(-0.1 * (L - x)); % Damping on the right side

    % Combine the left and right damping coefficients
    damping_coefficient = damping_coefficient_left .* damping_coefficient_right;

    % Initialize arrays for 3D plotting
    [X, T] = meshgrid(x, tspan(1):dt:tspan(2));
    U3D = zeros(N, num_frames + 1); % Updated dimensions
    U3D(:, 1) = v;

    % Precompute coefficients for the pseudo-spectral method
    k = fftshift((2 * pi / L) * (-N/2:N/2-1)); % Wavenumbers
    alpha = epsilon / 2;
    K = exp(-alpha * k.^2);

    % Create a VideoWriter object
    videoFileName = 'allen_cahn_animation.mp4';
    v = VideoWriter(videoFileName, 'MPEG-4');
    open(v);

    % Create a figure for 3D plot
    figure;
    ax = gca;
    ax.ZLim = [-1, 1]; % Adjust z-axis limits as needed

    % Initialize the 3D plot
    h = plot3(x, tspan(1) * ones(N, 1), U3D(:, 1));
    title('Allen-Cahn 1D Simulation');
    xlabel('x');
    ylabel('Time');
    zlabel('u');
    grid on;
    hold on;

    % Time-stepping loop using Crank-Nicolson
    for i = 1:num_frames
        U_current = U3D(:, i);

        % Crank-Nicolson time-stepping scheme with damping
        A = speye(N) + 0.5 * dt * epsilon * diag(k.^2) - dt * diag(damping_coefficient);
        B = speye(N) - 0.5 * dt * epsilon * diag(k.^2) + dt * diag(damping_coefficient);

        U_next = A \ (B * U_current);
        U3D(:, i + 1) = U_next;

        % Update the 3D plot
        set(h, 'XData', x, 'YData', (tspan(1) + i * dt) * ones(N, 1), 'ZData', U_next);
        title(['Time: ', num2str(tspan(1) + i * dt)]);
        
        drawnow; % Update the figure

        % Capture the frame and write it to the video
        frame = getframe(gcf);
        writeVideo(v, frame);

        pause(0.05); % Add a short pause for smoother animation
    end

    % Close the video file
    close(v);

    disp('Animation saved as "allen_cahn_animation.mp4".');
end
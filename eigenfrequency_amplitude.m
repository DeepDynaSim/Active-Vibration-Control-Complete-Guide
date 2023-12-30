% Adjusted Parameters 
E = 30e9; % Young's modulus (Pa)
I = 1e-8; % Second moment of inertia (m^4)
rho = 7800; % Density (kg/m^3)
S = 0.01; % Cross-sectional area (m^2)
L = 2; % Length of the beam (m)
gamma = 0.01; % Damping coefficient (1/s)
z0 = 0.01; % Amplitude of oscillation at z=0 (m)

% Calculating c
c = E*I / (rho*S);

% Frequency range
omega_range = linspace(0, 200, 10000); % 0-200 rad/s

% Vectorized Calculation of Amplitude
w_squared_over_c_squared = (omega_range.^2 / c^2);
q = sqrt(sqrt(w_squared_over_c_squared ./ (1 - 1i * gamma ./ omega_range)));
cos_qL = cos(q*L);
cosh_qL = cosh(q*L);
A_omega = abs((2*z0) ./ (2 * (1 + cosh_qL .* cos_qL)) .* (cos_qL + cosh_qL));

% Plotting Amplitude vs Frequency
figure;
plot(omega_range, A_omega);
xlabel('Frequency (rad/s)');
ylabel('Amplitude (m)');
title('Amplitude of the Free End of the Beam vs Frequency');

% Function for Transcendental Equation
Delta = @(q) 2 * (1 + cosh(q*L).*cos(q*L));

% Root Finding for Eigenfrequencies
qRange = linspace(0, 100, 10000);
qRoots = [];
for i = 1:length(qRange)-1
    if sign(Delta(qRange(i))) ~= sign(Delta(qRange(i+1)))
        root = fzero(Delta, [qRange(i), qRange(i+1)]);
        qRoots = [qRoots, root];
    end
end

% Eigenfrequencies Calculation
eigenfrequencies = sqrt(qRoots.^4 * c^2); % Vectorized Calculation

% Display Eigenfrequencies
disp('Eigenfrequencies (rad/s):');
disp(eigenfrequencies);

% Plotting for the First Four Resonant Frequencies
numResonantFreqs = min(length(eigenfrequencies), 4);
delta = 5; % Range around the resonant frequency
for i = 1:numResonantFreqs
    res_freq = eigenfrequencies(i);
    idx = omega_range > (res_freq - delta) & omega_range < (res_freq + delta);
    figure;
    plot(omega_range(idx), A_omega(idx));
    xlabel('Frequency (rad/s)');
    ylabel('Amplitude (m)');
    title(['Amplitude Near Resonant Frequency ', num2str(i)]);
end
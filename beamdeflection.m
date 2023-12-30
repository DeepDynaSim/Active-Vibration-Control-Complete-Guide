% Beam Deflection Analysis using Euler-Bernoulli Beam Theory

% Given Data
L = 6;          % Length of the beam in meters
E = 200e9;      % Modulus of Elasticity in Pa
I = 8e-6;       % Moment of Inertia in m^4
w = 5000;       % Uniformly Distributed Load in N/m
P = 10000;      % Point Load in N
w_max = 8000;   % Maximum Intensity of Triangular Load in N/m

% Calculate Reactions at Supports
R1 = (w * L / 2) + (P / 2) + (w_max * L / 2) / 2; % Reaction at left support
R2 = R1; % Reaction at right support (symmetric loading)

% Shear Force and Bending Moment Calculations
syms x C1 C2;

% Shear Force (V) and Bending Moment (M) as functions of x
V = R1 - w * x - (P / 2) * heaviside(x - L/2) - (w_max * x / L) * x / 2;
M = int(V, x);

% Deflection Analysis using Euler-Bernoulli Beam Equation
d2y = diff(M, x) / (E * I);
y_int = int(int(d2y, x), x) + C1 * x + C2;

% Apply boundary conditions to solve for integration constants
sol = solve(subs(y_int, x, 0) == 0, subs(y_int, x, L) == 0, C1, C2);

% Substitute back the constants
y = subs(y_int, [C1, C2], [sol.C1, sol.C2]);

% Define a function handle for numerical evaluation
y_func = matlabFunction(y);

% Calculate Deflection at Specific Points
deflection_midpoint = double(y_func(L/2));
deflection_2m = double(y_func(2));
deflection_right_support = double(y_func(L));

% Find Maximum Deflection
x_max_deflection = fminbnd(@(x) -y_func(x), 0, L);
max_deflection = -y_func(x_max_deflection);

% Display Results
fprintf('Reactions at Supports: R1 = %.2f N, R2 = %.2f N\n', R1, R2);
fprintf('Deflection at Midpoint: %.4e m\n', deflection_midpoint);
fprintf('Deflection at 2m from Left Support: %.4e m\n', deflection_2m);
fprintf('Deflection at Right Support: %.4e m\n', deflection_right_support);
fprintf('Maximum Deflection: %.4e m at x = %.2f m\n', max_deflection, x_max_deflection);

% Plotting Beam Analysis
figure;
subplot(4,1,1);
fplot(V, [0, L], 'b', 'LineWidth', 2);
title('Shear Force Diagram');
ylabel('V (N)');

subplot(4,1,2);
fplot(M, [0, L], 'r', 'LineWidth', 2);
title('Bending Moment Diagram');
ylabel('M (N.m)');

subplot(4,1,3);
fplot(y_func, [0, L], 'g', 'LineWidth', 2);
title('Deflection Curve');
xlabel('x (m)');
ylabel('Deflection (m)');

% Visual Representation of Beam and Loading Conditions
figure(2)
hold on;
grid on;

% Draw the Beam
plot([0, L], [0, 0], 'k', 'LineWidth', 2);

% Draw Supports
fill([0, 0.2, -0.2, 0], [0, -0.2, -0.2, 0], 'b'); % Left Support
fill([L, L+0.2, L-0.2, L], [0, -0.2, -0.2, 0], 'b'); % Right Support

% Draw Loads
% UDL
for i = 0:0.2:L
    plot([i, i], [0, -0.3], 'g'); % UDL lines
end
% Point Load
plot([L/2, L/2], [0, -0.5], 'r', 'LineWidth', 2); % Point load line
% Triangular Load
fill([0, L, L], [-0.6, -0.6, 0], 'm'); % Triangular load

% Annotations
text(L/2, -0.15, 'UDL', 'HorizontalAlignment', 'center', 'Color', 'g');
text(L/2, -0.4, 'P', 'HorizontalAlignment', 'center', 'Color', 'r');
text(0, -0.65, 'Triangular Load', 'Color', 'm');

% Set Axis Limits and Labels
axis([-1 L+1 -1 0.5]);
xlabel('Length (m)');
ylabel('Load (N or N/m)');
title('Beam with Supports and Loading Conditions');

hold off;
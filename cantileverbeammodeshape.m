% Cantilever Beam Vibration Problem
format long e;

% Define symbolic variable
syms x real;

% Parameters for Cantilever Beam Vibration Problem
E = 180e9;      % Elastic Modulus (N/m^2)
rho = 7500;     % Mass density (kg/m^3)
L = 3.9;        % Length (m)
A = 0.001;      % Cross-sectional area (m^2) - assuming constant for simplicity
I = A^2 / 12;   % Area moment of inertia (m^4) - assuming a rectangular cross section

% Define beta based on the problem
% Replace this with the specific expression
beta = sqrt(sqrt(E*I/(rho*A))/L); % Example expression

% Define mode shapes based on the problem
% Replace with the specific mode shape expressions
modeShapes = [sinh(beta*x) - sin(beta*x) - cosh(beta*x) + cos(beta*x), ...
              sinh(2*beta*x) - sin(2*beta*x) - cosh(2*beta*x) + cos(2*beta*x), ...
              sinh(3*beta*x) - sin(3*beta*x) - cosh(3*beta*x) + cos(3*beta*x)];
disp('Solving Cantilever Beam Vibration Problem');

% Preallocating matrices for Mass (M) and Stiffness (K)
M = sym(zeros(3));
K = sym(zeros(3));

% Calculate Mass and Stiffness matrices
for i = 1:3
    for j = 1:i
        Mint = rho * A * modeShapes(i) * modeShapes(j);
        Kint = E * A * diff(modeShapes(i), x, 2) * diff(modeShapes(j), x, 2);
        M(i, j) = int(Mint, x, 0, L);
        K(i, j) = int(Kint, x, 0, L);
        % Symmetric matrices
        M(j, i) = M(i, j);
        K(j, i) = K(i, j);
    end
end

% Converting to double for numerical calculations
K_num = double(K);
M_num = double(M);

% Calculating the system matrix
C = inv(M_num) * K_num;

% Eigenvectors and eigenvalues
[V, D] = eig(C);

% Natural frequencies
w = sqrt(diag(D));

% Normalizing eigenvectors (Modal Matrix)
P = V' * M * V;
for j = 1:3
    for i = 1:3
        P(i, j) = V(i, j) / sqrt(P(j, j));
    end
end

% Display results
disp('Natural frequencies in rad/s');
disp(vpa(w));
disp('Modal matrix');
disp(vpa(P));

% Mode shape visualization
xx = linspace(0, L, 100);
v = zeros(3, length(xx));
for k = 1:length(xx)
    for i = 1:3
        for j = 1:3
            v(i, k) = v(i, k) + P(j, i) * subs(modeShapes(j), x, xx(k));
        end
    end
end

% Plot mode shapes
figure;
plot(xx, v(1, :), '-', xx, v(2, :), '-.', xx, v(3, :), '--');
xlabel('x (m)');
ylabel('w(x) (m)');
legend('1st mode', '2nd mode', '3rd mode');
title('Mode Shapes of the Cantilever Beam');
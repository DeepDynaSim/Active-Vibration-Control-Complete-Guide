%The provided MATLAB code is a well-structured script 
% for Finite Difference Method (FDM) analysis of a beam under uniform load.
% Finite Difference Method Analysis of a Beam Under Uniform Load
% Finite Difference Method Analysis of a Beam Under Uniform Load

function beamdeflection_fdm(L, I, E, w0, ne)
    % Initialization
    EI = E * I;    % Flexural rigidity
    nx = ne + 1;   % Number of nodes
    dx = L / ne;   % Size of elements
    X = 0:dx:L;    % Node positions along the beam

    % Load vector
    f = -w0 * dx^4 / EI * ones(nx, 1);

    % Stiffness matrix (sparse for efficiency)
    Kstiff = sparse(nx, nx);

    % General Differential Equation (GDE) Implementation
    Kstiff = buildStiffnessMatrix(Kstiff, nx);

    % Apply Boundary Conditions
    [Kstiff, f] = applyBoundaryConditions(Kstiff, f, nx);

    % Solve for deflection
    y = Kstiff \ f;

    % Calculate moment, shear, and pressure
    [M, V, P] = computeDifferentialValues(y, EI, dx, nx);

    % Display maximum values
    displayMaxValues(y, M, V, P);

    % Plot results
    plotResults(X, y, M, V, P);
end

function Kstiff = buildStiffnessMatrix(Kstiff, nx)
    diagonal = [1 -4 6 -4 1];
    for i = 3:nx-2
        Kstiff(i, i-2:i+2) = diagonal;
    end
end

function [Kstiff, f] = applyBoundaryConditions(Kstiff, f, nx)
    Kstiff(1,1) = 1;
    f(1) = 0;
    Kstiff(2,1:3) = [1 -2 1];
    Kstiff(nx, nx) = 1;
    f(nx) = 0;
    Kstiff(nx-1, nx-2:nx) = [1 -2 1];
end

function [M, V, P] = computeDifferentialValues(y, EI, dx, nx)
    M = zeros(nx, 1);
    V = zeros(nx, 1);
    P = zeros(nx, 1);

    for i = 2:nx-1
        M(i) = computeMoment(y, i, EI, dx);
        if i > 2 && i < nx-1
            V(i) = computeShear(y, i, EI, dx);
            P(i) = computePressure(y, i, EI, dx);
        end
    end
end

function M = computeMoment(y, i, EI, dx)
    M = EI * (y(i-1) - 2*y(i) + y(i+1)) / dx^2;
end

function V = computeShear(y, i, EI, dx)
    V = EI * (-y(i-2) + 3*y(i-1) - 3*y(i) + y(i+1)) / dx^3;
end

function P = computePressure(y, i, EI, dx)
    P = EI * (y(i-2) - 4*y(i-1) + 6*y(i) - 4*y(i+1) + y(i+2)) / dx^4;
end

function displayMaxValues(y, M, V, P)
    fprintf('Max displacement = %f \n', max(abs(y)));
    fprintf('Max moment = %f \n', max(abs(M)));
    fprintf('Max shear = %f \n', max(abs(V)));
    fprintf('Max pressure = %f \n', max(abs(P)));
end

function plotResults(X, y, M, V, P)
    fontSize = 20;
    lineWidth = 2;

    subplot(4, 1, 1);
    plotBeam(X, y, 'Deflection', 'b', fontSize, lineWidth);
    subplot(4, 1, 2);
    plotBeam(X, M, 'Moment', 'r', fontSize, lineWidth);
    subplot(4, 1, 3);
    plotBeam(X, V, 'Shear', 'g', fontSize, lineWidth);
    subplot(4, 1, 4);
    plotBeam(X, P, 'Pressure', 'm', fontSize, lineWidth);
end

function plotBeam(X, Y, label, color, fontSize, lineWidth)
    stem(X, Y, color, 'LineWidth', lineWidth);
    xlabel('Beam Length', 'FontSize', fontSize);
    ylabel(label, 'FontSize', fontSize);
    set(gca, 'FontSize', fontSize);
    grid on;
    title(sprintf('FDM Analysis for Beam with Uniform Load - %s', label), 'FontSize', fontSize);
end

%Example of Realistic Parameters:
%Length of the beam, L: 10 meters
%Moment of inertia, I: 0.0004 
%Elastic modulus, E: 200 GPa (200 x 10^9 N/m^2)
%Uniformly distributed load, w0: 5000 N/m
%Number of elements, ne: 50
%beamdeflection_fdm(10, 0.0004, 200e9, 5000, 50);
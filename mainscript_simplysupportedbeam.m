% Number of elements
m = 10; 

% Load in Newtons
% Assuming a uniformly distributed load (UDL)
P = -1000; % Negative sign indicates downward load

% Beam length in meters
L = 10; % For example, a 10-meter long beam

% Modulus of Elasticity in Pascals (e.g., for steel approximately 200 GPa)
E = 200e9; % 200 GPa converted to Pascals

% Moment of Inertia in meters^4
% This depends on the cross-sectional shape and size of the beam
% For example, for a rectangular cross-section (b x h), I = (b*h^3)/12
I = 1e-4; % Example value, adjust based on actual cross-section

% Call the function with the defined parameters
[stiffness, force, displacements, reactions] = formStiffness_simplysupportedbeam(m, P, L, E, I);
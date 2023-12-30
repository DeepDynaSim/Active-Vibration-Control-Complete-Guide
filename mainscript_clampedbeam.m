% Main script
E = 210e9; % Modulus of elasticity for steel (in Pascals)
I = 1e-6; % Moment of inertia (in m^4)
m = 10; % Number of elements

% Plot displacements under various loads
loads = -25:-25:-100; % Loads from 25 N to 100 N
figure;
hold on;
for P = loads
    displacements = formStiffness_clampedbeam(m, P, E, I);
    plot(displacements(1:2:end), '-o', 'DisplayName', sprintf('Load = %d N', P));
end
title('Displacements of a Clamped Beam under Various Loads');
xlabel('Node Number');
ylabel('Displacement (m)');
legend show;
hold off;

% Plot displacements for different element counts under a specific load
elementCounts = 10:5:30; % Element counts from 10 to 30
P = -1000; % Specific load
figure;
hold on;
for m = elementCounts
    displacements = formStiffness_clampedbeam(m, P, E, I);
    plot(displacements(1:2:end), '-o', 'DisplayName', sprintf('%d Elements', m));
end
title('Displacements of a Clamped Beam with Different Element Counts');
xlabel('Node Number');
ylabel('Displacement (m)');
legend show;
hold off;
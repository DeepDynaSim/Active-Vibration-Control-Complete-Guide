function displacements = formStiffness_clampedbeam(m, P, E, I)
    % Initialize parameters
    L = 1; % Length of the beam
    nodeCoordinates = linspace(0, L, m+1)'; 
    numberNodes = numel(nodeCoordinates);
    EI = E * I;
    GDof = 2 * numberNodes; % Total number of degrees of freedom

    % Initialize matrices and vectors
    stiffness = zeros(GDof);
    force = zeros(GDof, 1);

    % Element connectivity and assembly
    for e = 1:m
        elementNodes = [e, e + 1];
        elementDof = [2*elementNodes(1)-1, 2*elementNodes(1), 2*elementNodes(2)-1, 2*elementNodes(2)];
        LElem = L / m; % Length of each element
        f1 = P * LElem / 2 * [1, LElem/6, 1, -LElem/6]';
        force(elementDof) = force(elementDof) + f1;
        k1 = EI / LElem^3 * [12, 6*LElem, -12, 6*LElem; 6*LElem, 4*LElem^2, -6*LElem, 2*LElem^2; -12, -6*LElem, 12, -6*LElem; 6*LElem, 2*LElem^2, -6*LElem, 4*LElem^2];
        stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + k1;
    end

    % Apply boundary conditions
    prescribedDof = [1; 2; 2*m+1; 2*m+2];
    activeDof = setdiff((1:GDof)', prescribedDof);

    % Solve for displacements
    displacements = zeros(GDof, 1);
    displacements(activeDof) = stiffness(activeDof, activeDof) \ force(activeDof);
end
function [stiffness, force, displacements, reactions] = formStiffness_simplysupportedbeam(m, P, L, E, I)
    % Parameters
    nodeCoordinates = linspace(0, L, m+1)'; % Beam length in meters
    numberNodes = size(nodeCoordinates, 1);
    EI = E * I; % Modulus of elasticity (Pa) * Moment of inertia (m^4)
    GDof = 2 * numberNodes; % Total number of degrees of freedom

    % Initialize matrices and vectors
    force = zeros(GDof, 1);
    stiffness = zeros(GDof);
    displacements = zeros(GDof, 1);

    % Element connectivity
    elementNodes = zeros(m, 2);
    for i = 1:m
        elementNodes(i, 1) = i;
        elementNodes(i, 2) = i + 1;
    end

    % Assembly of stiffness matrix and force vector
    for e = 1:m
        indice = elementNodes(e, :);
        elementDof = [2*(indice(1)-1)+1, 2*indice(1), 2*(indice(2)-1)+1, 2*indice(2)];
        LElem = nodeCoordinates(indice(2)) - nodeCoordinates(indice(1));
        
        % Uniformly Distributed Load (UDL) assumption
        f1 = P * LElem * [1/2, LElem/6, 1/2, -LElem/6]';
        force(elementDof) = force(elementDof) + f1;
        
        % Stiffness matrix for each element
        k1 = EI/(LElem)^3 * [12, 6*LElem, -12, 6*LElem; 6*LElem, 4*LElem^2, -6*LElem, 2*LElem^2; -12, -6*LElem, 12, -6*LElem; 6*LElem, 2*LElem^2, -6*LElem, 4*LElem^2];
        stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + k1;
    end

    % Apply boundary conditions for simply supported beam
    fixedNodeU = [1, 2*m+1]'; % Vertical displacement constraints at both ends
    fixedNodeV = []; % No rotational fixity constraints
    prescribedDof = [fixedNodeU; fixedNodeV];
    activeDof = setdiff([1:GDof]', prescribedDof);

    % Solve for displacements
    U = zeros(length(activeDof), 1);
    if ~isempty(stiffness(activeDof, activeDof))
        U = stiffness(activeDof, activeDof) \ force(activeDof);
    end
    displacements(activeDof) = U;

    % Calculate reactions
    F = stiffness * displacements;
    reactions = F(prescribedDof);

    % Display results
    disp('Displacements');
    plot(nodeCoordinates, displacements(1:2:GDof), '-o');
    xlabel('Position along Beam');
    ylabel('Vertical Displacement');
    title('Displacement of Simply Supported Beam under Load');
end
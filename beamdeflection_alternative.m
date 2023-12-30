function [endDeflection, endSlope, DeflectionVector] = beamdeflection_alternative(loadForce, E, beamLength, beamWidth, type)
    % BEAMDEFLECTION Calculate deflection of a cantilever beam under load
    %
    % INPUTS:
    % loadForce     - Concentrated load [N]
    % E             - Young's modulus of the material [N/m^2]
    % beamLength    - Total length of the beam [m]
    % beamWidth     - Width of the beam [m]
    % type          - Cross-section shape; 'rectangle', 'circle', 'triangle', 'hexagon'
    %
    % OUTPUTS:
    % endDeflection - Deflection of the end-tip of the beam [m]
    % endSlope      - Angle with the horizontal at the end-tip of the beam [-]
    % DeflectionVector - Vector of deflection distances along the beam [m]

    % Check for sufficient input arguments
    if nargin < 4
        error('Insufficient input arguments. Please provide loadForce, E, beamLength, and beamWidth.');
    end

    % Default to 'rectangle' if no type specified
    if nargin < 5
        type = 'rectangle';
        warning('beamDeflection:defaultType', 'No beam type specified. Using "rectangle"');
    end

    % Validate inputs
    validateattributes(loadForce, {'numeric'}, {'scalar', 'positive'});
    validateattributes(E, {'numeric'}, {'scalar', 'positive'});
    validateattributes(beamLength, {'numeric'}, {'scalar', 'positive'});
    validateattributes(beamWidth, {'numeric'}, {'scalar', 'positive'});
    type = validatestring(type, {'rectangle', 'circle', 'triangle', 'hexagon'});

    % Calculate moment of inertia, Ixx, based on cross-sectional shape
    Ixx = calcMomentOfInertia(beamWidth, type);

    % Compute beam deformation
    [endDeflection, endSlope, DeflectionVector] = calcDeformation(loadForce, E, Ixx, beamLength);

    % Plot results if no output arguments
    if nargout == 0
        plotBeamDeflection(beamLength, DeflectionVector, type, endDeflection);
    end
end

function Ixx = calcMomentOfInertia(beamWidth, type)
    % Calculate moment of inertia for different cross-sectional shapes
    switch type
        case 'rectangle'
            Ixx = (beamWidth^4) / 12;
        case 'circle'
            Ixx = (pi * beamWidth^4) / 64;
        case 'triangle'
            Ixx = (beamWidth^4) / 36;
        case 'hexagon'
            Ixx = (5 * sqrt(3) * (0.5 * beamWidth)^4) / 16;
    end
end

function [endDeflection, endSlope, DeflectionVector] = calcDeformation(loadForce, E, Ixx, beamLength)
    % Calculate beam deformation parameters
    endDeflection = (loadForce * beamLength^3) / (3 * E * Ixx);
    endSlope = (loadForce * beamLength^2) / (2 * E * Ixx);
    beam = linspace(0, beamLength, 100);
    DeflectionVector = (loadForce * beam.^2) / (6 * E * Ixx) .* (3 * beamLength - beam);
end

function plotBeamDeflection(beamLength, DeflectionVector, type, endDeflection)
    % Plot the beam deflection
    beam = linspace(0, beamLength, 100);
    figure
    plot(beam, zeros(size(beam)), 'k--')
    hold on
    plot(beam, -DeflectionVector, 'LineWidth', 2)
    ylim([-3 * endDeflection, 3 * endDeflection])
    xlim([0, 1.1 * beam(end)])
    title(['Cross-section: ', type])
    xlabel('Distance from fixed support [m]')
    ylabel('Deflection [m]')
end

%[endDeflection, endSlope, DeflectionVector] = beamdeflection_alternative(1500, 210e9, 3, 0.15, 'circle');


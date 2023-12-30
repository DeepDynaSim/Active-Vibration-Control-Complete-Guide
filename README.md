# Active-Vibration-Control-Complete-Guide

A Comprehensive Guide to Active Control Methodology for Smart Beams
This repository houses the groundbreaking research and development of a methodology for active vibration control in Euler-Bernoulli-type smart beams. Our project merges intricate numerical analysis with robust-based boundary control synthesis. The primary aim is to significantly attenuate unwanted vibrations in engineering structures by harnessing the capabilities of piezoelectric-based smart beam technology.
Key Features:
Numerical Analysis: We focus into the numerical analysis of smart beams, focusing on the accurate simulation of the Euler-Bernoulli smart beam model. This analysis is essential for effectively modeling different beam configurations under a variety of vibration conditions. It's the cornerstone for determining optimal control parameters and appropriate piezoelectric actuation voltages.
Robust Boundary Control Synthesis: The crux of this repository is the introduction of a robust-based boundary control approach. Our system design synergizes boundary control techniques with robust algorithms, paving the way for an efficient control mechanism. This hybrid strategy enables the system to learn and evolve its control precision autonomously, dynamically responding to environmental fluctuations.
Validation through Realistic Simulations: The methodology's effectiveness is corroborated through extensive, realistic simulations. These simulations affirm the capacity of our approach to suppress vibrations, with results closely mirroring our numerical and analytical projections. This not only validates the robustness of our integrated approach but also its practical applicability in real-world scenarios.
Applications:
The methodology and findings presented in this repository are pivotal for professionals and researchers in the fields of structural engineering, robotics, aerospace, and any domain where active vibration control is crucial. Our approach sets a new benchmark in smart beam technology and offers a sophisticated toolkit for tackling vibration-related challenges in engineering structures.
Contributions and Collaborations:
We welcome contributions from fellow researchers, engineers, and enthusiasts. Whether it's refining algorithms, enhancing simulation models, or discussing potential applications, your insights and expertise can help push the boundaries of this innovative technology.

beamdeflection.m 

Description:

This repository provides a MATLAB script for analyzing beam deflection using the Euler-Bernoulli beam theory. The script is a comprehensive tool for engineering students and professionals to understand and compute the deflection, shear force, and bending moment in beams subjected to various load conditions. The code is well-commented, making it easy to follow and adapt for specific beam analysis tasks.

Step-by-Step Explanation:

Initial Setup:

Data Definition: The script starts by defining the beam's physical properties and load conditions, including length (L), modulus of elasticity (E), moment of inertia (I), uniformly distributed load (w), point load (P), and maximum intensity of a triangular load (w_max).
Reaction Forces Calculation:

Support Reactions: It calculates the reaction forces at the beam's supports (R1, R2), considering the symmetry of the loading.
Shear Force and Bending Moment Analysis:

Symbolic Variables: Defines symbolic variables x, C1, and C2 for integration and differentiation.
Shear Force (V): Computes the shear force as a function of x, incorporating various loads (uniform, point, and triangular).
Bending Moment (M): Calculates the bending moment by integrating the shear force over x.
Deflection Analysis:

Euler-Bernoulli Equation: Applies the Euler-Bernoulli beam equation to find the beam's deflection. This involves differentiating the bending moment and integrating twice.
Boundary Conditions: Solves for the integration constants C1 and C2 using the boundary conditions (deflection is zero at both ends of the beam).
Numerical Evaluation and Maximum Deflection:

Function Handle: Creates a MATLAB function handle for deflection (y_func) to facilitate numerical evaluation.
Specific Deflections: Calculates deflections at the midpoint, 2 meters from the left support, and at the right support.
Maximum Deflection: Determines the point of maximum deflection on the beam and its value.
Results Display:

Printing Outputs: Displays the calculated reaction forces, deflections at specific points, and maximum deflection.
Graphical Representation:

Plotting: Generates plots for the shear force diagram, bending moment diagram, and deflection curve.
Visualizing Beam and Loads: Draw a schematic of the beam with supports and applied loads, enhancing the understanding of the load distribution and its impact.

beamdeflection_alternative.m

Overview:
Function Purpose: To calculate and potentially visualize the deflection of a cantilever beam under various conditions.
Inputs: Load force, Young's modulus, beam length, beam width, and cross-section type.
Outputs: End-tip deflection, end-tip slope angle, and a vector of deflection distances along the beam.
Function Breakdown:
Input Validation and Defaults:

Checks if sufficient arguments are provided. If not, it throws an error.
Sets the default cross-section type to 'rectangle' if not specified.
Validates the input types and constraints using validate attributes and validate string.
Moment of Inertia Calculation (calcMomentOfInertia):

This function calculates the moment of inertia (Ixx) based on the beam's cross-sectional shape. It handles four shapes: rectangle, circle, triangle, and hexagon.
The moment of inertia is essential for determining the beam's resistance to bending.
Deformation Calculation (calcDeformation):

This function calculates the end-tip deflection, end-tip slope, and a vector of deflections along the beam length.
These calculations use standard formulas from beam theory, incorporating the moment of inertia, beam properties, and load force.
Result Visualization (plotBeamDeflection):

If no output arguments are specified, the function plots the beam deflection.
The plot includes a line representing the beam and its deflection curve, with labels and a title indicating the cross-section type.
Example Usage:
The commented-out line at the end shows an example call to the function with specific parameters (load force of 1500 N, Young's modulus of 210e9 N/m^2, beam length of 3 m, beam width of 0.15 m, and a circular cross-section).
This example would calculate and potentially plot the deflection of a circular cantilever beam under the specified conditions.
Educational and Practical Significance:
For Students: This function is an excellent tool for learning how different parameters affect beam deflection. It provides a hands-on approach to understanding cantilever beam mechanics.
For Experts: It serves as a quick and efficient way to estimate cantilever beam deflections for different cross-sections, which is useful in the preliminary design and analysis phases.
Extensibility:
The function is designed to be easily extendable. Additional cross-sectional shapes or more complex loading conditions can be incorporated with minimal modifications.

beamdeflection_fdm.m

Purpose: Main function to perform FDM analysis of a uniformly loaded beam.
Inputs: Beam length L, moment of inertia I, elastic modulus E, uniformly distributed load w0, and number of elements ne.
Steps in beamdeflection_fdm:
Initialization: Computes flexural rigidity (EI), the number of nodes (nx), the size of elements (dx), and the node positions (X).
Load Vector: Creates a load vector f representing the distributed load along the beam.
Stiffness Matrix: Initializes and builds a sparse stiffness matrix Kstiff using buildStiffnessMatrix.
Boundary Conditions: Applies boundary conditions to the stiffness matrix and load vector.
Solve for Deflection: Solves the linear system to find the deflection at each node.
Compute Differential Values: Calculates moment (M), shear (V), and pressure (P) along the beam.
Display Maximum Values: Prints maximum values of displacement, moment, shear, and pressure.
Plot Results: Visualizes the results for deflection, moment, shear, and pressure.
Supporting Functions:
buildStiffnessMatrix: Constructs the stiffness matrix for the FDM analysis.
apply boundary conditions: Modifies the stiffness matrix and load vector to incorporate boundary conditions.
computeDifferentialValues: Calculates moment, shear, and pressure at each node.
computeMoment, computeShear, computePressure: Helper functions to compute moment, shear force, and pressure, respectively, at a specific node.
displayMaxValues: Displays the maximum values of displacement, moment, shear, and pressure.
plotResults: Plots the deflection, moment, shear, and pressure along the beam.
plotBeam: A utility function to plot a specific parameter along the beam.
Example Usage:
The code is designed to be used with parameters like a 10-meter beam (L), a moment of inertia of 0.0004 (I), an elastic modulus of 200 GPa (E), a uniformly distributed load of 5000 N/m (w0), and divided into 50 elements (ne).
This realistic setup would allow users to analyze a beam's response under specified loading conditions, providing valuable insights into its structural behavior.


The MATLAB function formStiffness_clampedbeam is designed to calculate the displacements of a clamped beam under a uniform load using the finite element method (FEM). Let's break down the code into its key components:

Function Overview:
Purpose: To compute the displacements of a clamped beam subjected to a uniform load.
Inputs: Number of elements m, load per unit length P, Young's modulus E, and moment of inertia I.
Output: A vector of displacements at each degree of freedom in the beam.
Steps in formStiffness_clampedbeam:
Initialization:

Sets the length of the beam (L) to 1 meter.
Generates node coordinates evenly spaced along the beam.
Calculates the total number of degrees of freedom (GDof), considering two degrees of freedom per node (for bending in beams, typically vertical displacement and rotation).
Matrix and Vector Initialization:

Initializes the stiffness matrix (stiffness) and force vector (force) with zeros. The size of these matrices is based on the total degrees of freedom.
Element Connectivity and Assembly:

Loops over each element (m) to calculate and assemble the local stiffness matrix (k1) and force vector (f1) for each element.
For each element, it identifies the local degrees of freedom (elementDof) and updates the global stiffness matrix and force vector accordingly.
The local stiffness matrix is derived from beam theory (Euler-Bernoulli beam theory) and accounts for the bending stiffness of the beam segments.
Boundary Conditions:

Identifies the degrees of freedom that are constrained due to clamped conditions at both ends of the beam (prescribedDof).
Determines the active degrees of freedom (activeDof) where the displacements are unknown and need to be solved.
Solving for Displacements:

Initializes the displacement vector (displacements) with zeros.
Solves the system of linear equations for the active degrees of freedom to find the displacements. This is done by restricting the stiffness matrix and force vector to the active degrees of freedom and using matrix division.

The MATLAB function formStiffness_simplysupportedbeam calculates the displacements and reactions of a simply supported beam under a uniformly distributed load using the Finite Element Method (FEM). Here's a detailed breakdown of its functionality:

Function Overview:
Purpose: To compute displacements and reactions for a simply supported beam under a uniform load.
Inputs: Number of elements m, load per unit length P, beam length L, Young's modulus E, and moment of inertia I.
Outputs: Stiffness matrix (stiffness), force vector (force), displacement vector (displacements), and reaction forces (reactions).
Steps in formStiffness_simplysupportedbeam:
Parameter Initialization:

Sets up node coordinates along the beam based on the specified length L and number of elements m.
Calculates the total degrees of freedom (GDof), accounting for both vertical displacement and rotation at each node.
Matrix and Vector Initialization:

Initializes the force vector, stiffness matrix, and displacement vector with zero values.
Element Connectivity:

Establishes the connectivity of the elements, identifying which nodes are connected by each element.
Stiffness Matrix and Force Vector Assembly:

Loops over each element, computing and assembling the local stiffness matrix and force vector for each element.
The local stiffness matrix (k1) reflects the bending behavior of the beam elements.
The force vector (f1) is calculated assuming a uniformly distributed load.
Boundary Conditions Application:

Applies boundary conditions specific to a simply supported beam, i.e., vertical displacement constraints at both ends of the beam.
Identifies active degrees of freedom where displacements are unknown.
Solving for Displacements:

Solves the linear system for displacements at the active degrees of freedom.
Updates the displacement vector with the computed values.
Reaction Forces Calculation:

Multiplies the stiffness matrix by the displacement vector to calculate the forces, and extracts the reaction forces at the prescribed degrees of freedom.
Results Display and Visualization:

Displays the calculated displacements graphically, showing vertical displacement along the beam.

cantileverbeammodeshape.m

The MATLAB script provided performs a vibration analysis of a cantilever beam, a fundamental problem in structural dynamics and vibration engineering. The script uses a combination of symbolic and numerical computation to find natural frequencies and mode shapes of the beam. Let's examine the script in detail:

Script Overview:
Purpose: To analyze the vibration characteristics (natural frequencies and mode shapes) of a cantilever beam.
Key Parameters: Elastic modulus E, mass density rho, beam length L, cross-sectional area A, and area moment of inertia I.
Steps in the Script:
Parameter Initialization:

Defines the physical properties of the beam, such as elastic modulus, mass density, length, cross-sectional area, and area moment of inertia.
Beta Calculation:

Calculates beta, a parameter used in the definition of mode shapes, based on the beam's properties.
Mode Shapes Definition:

Defines the first three mode shapes of the cantilever beam as symbolic expressions in terms of x.
Mass and Stiffness Matrices Computation:

Initializes and calculates the mass (M) and stiffness (K) matrices by integrating the products of mode shapes and their derivatives. These matrices are essential in vibration analysis.
Numerical Conversion and System Matrix Calculation:

Converts the mass and stiffness matrices to numerical values (M_num, K_num).
Computes the system matrix (C) by multiplying the inverse of the mass matrix with the stiffness matrix.
Eigenvalue Problem and Natural Frequencies:

Solves the eigenvalue problem for the system matrix to find the eigenvectors (V) and eigenvalues (D).
Calculates natural frequencies (w) as the square root of the eigenvalues.
Modal Matrix Normalization:

Normalizes the eigenvectors to form a modal matrix (P), ensuring that the mode shapes are orthogonal with respect to the mass matrix.
Results Display:

Outputs the natural frequencies and the modal matrix.
Mode Shape Visualization:

Evaluates the mode shapes over the length of the beam (xx) and plots them.
Plotting Mode Shapes:

Displays the first three mode shapes graphically, illustrating how the beam would vibrate in each mode.

The MATLAB function beam_vibration_ode45_armature simulates the dynamic response of a cantilever beam subject to an impulse load at its tip using numerical methods. The script is structured to model the physical behavior of the beam, solve the governing differential equations, and visualize the beam's response over time. Let's break down the script into its key components:

Script Overview:
Purpose: To analyze and visualize the vibration of a cantilever beam under an impulse load using numerical integration.
Key Elements: Elastic modulus E, mass density rho, beam dimensions L, b, h, and the impulse load characteristics.
Steps in the Script:
Parameter Initialization:

Sets the physical properties of the beam, including its material properties (elastic modulus and mass density), geometric dimensions (length, width, depth), area moment of inertia, and cross-sectional area.
Matrix Formulation:

Constructs the mass (M), stiffness (K), and damping (C) matrices. The stiffness matrix is adjusted for a cantilever beam, and Rayleigh damping is used to model energy dissipation.
Time Span and Initial Conditions:

Specifies the time span for the simulation and sets the initial conditions (zero displacement and velocity).
ODE Solver (ode45):

Uses MATLAB's ode45 function to numerically solve the ordinary differential equations (ODEs) governing the beam's vibration. The ODEs are defined in the beam_ode_impulse function.
Fixed-End Boundary Conditions:

Ensures that the displacement at the fixed end remains zero throughout the simulation.
Deflection Visualization:

Plots the deflection of the beam tip over time to visualize the vibration response.
Animation Creation:

Generates an animation of the beam's response to the impulse load. The animation is saved as an AVI file.
beam_ode_impulse Function:
Purpose: To define the system of ODEs for the beam's dynamic response.
Key Operation: Incorporates an impulse load at the beam tip, represented by a high magnitude force applied for a very short duration.
plot_beam Function:
Purpose: To visualize the displacement of the beam at a specific time instance.
Visualization: Scales the displacement for better visibility and plots the deformed shape of the beam.

forcedbeamvibration.m

The MATLAB script provided performs a numerical simulation of a cantilever beam undergoing vibrations. The beam's dynamic response to a harmonic driving force at one end is calculated using a discretized spatial domain and time-stepping method. Let's break down the key components of the script:

Script Overview:
Purpose: To simulate and analyze the vibrations of a cantilever beam under a harmonic driving force.
Key Parameters: Elastic modulus E, second moment of inertia I, density rho, cross-sectional area S, beam length L, damping coefficient gamma, amplitude of oscillation z0, and frequency of oscillation w.
Steps in the Script:
Parameter Initialization:

Sets the material and geometric properties of the beam, and the properties of the oscillation.
Discretization:

Defines the spatial (dz) and temporal (dt) discretization for the simulation.
Initial Conditions:

Sets initial deflection Fz and initial velocity Gz of the beam.
Time-Stepping Simulation:

Initializes the displacement matrix u and applies initial conditions.
Uses a combination of central difference in space and forward Euler in time for the time-stepping process.
Applies boundary conditions at each time step: a driven end with a harmonic force and a free end.
Plotting Results:

Visualizes the vibration of the beam over time using a mesh plot.
Calculation of c:

Calculates a constant c used in the eigenfrequency analysis.
Eigenfrequency Analysis:

Solves a transcendental equation to find roots (qRoots) using the function Delta.
Calculates the eigenfrequencies of the beam based on the roots found.
Display Eigenfrequencies:

Outputs the calculated eigenfrequencies.

eigenfrequency_amplitude.m


This MATLAB script is designed to perform a frequency response analysis of a cantilever beam. It involves calculating the response amplitude at different frequencies and identifying the eigenfrequencies of the beam. Let's break down the script into its key components:

Script Overview:
Purpose: To analyze the amplitude response of a cantilever beam to a range of frequencies and to determine the beam's eigenfrequencies.
Key Parameters: Elastic modulus E, second moment of inertia I, density rho, cross-sectional area S, beam length L, damping coefficient gamma, and amplitude z0 at z = 0.
Steps in the Script:
Parameter Initialization:

Sets the material and geometric properties of the beam, along with the damping coefficient and the amplitude of oscillation at one end.
Calculating c:

Computes the constant c, a parameter combining the material and geometric properties of the beam.
Frequency Response Calculation:

Creates a frequency range (omega_range) and calculates the beam's response amplitude (A_omega) at each frequency. This involves computing a complex wave number (q) and using it to calculate the amplitude.
Plotting Amplitude vs. Frequency:

Visualizes the amplitude response of the beam across the specified frequency range.
Eigenfrequency Analysis:

Uses a transcendental equation (Delta) to find the roots (qRoots), which correspond to the eigenfrequencies of the beam.
Vectorized calculation of the eigenfrequencies based on the roots found.
Display Eigenfrequencies:

Outputs the calculated eigenfrequencies.
Plotting Resonant Frequencies:

Focuses on the first four resonant frequencies and plots the amplitude response near these frequencies to illustrate resonance phenomena.

bernoullieuler_beamvibration_FDM.m

The provided MATLAB script is designed for simulating the vibration of a cantilever beam using the Finite Difference Method (FDM). The script includes parameters for the physical properties of the beam, applies boundary conditions, and visualizes the displacement of the beam over time. Here's a detailed breakdown of its key components:

Script Overview:
Purpose: To simulate and visualize the dynamic behavior of a cantilever beam under various loading conditions and boundary disturbances.
Key Parameters: Includes physical properties of the beam like length L, width wd, thickness t, elastic modulus E, density rho, and others.
Steps in the Script:
Parameter Initialization:

Sets the physical and geometric properties of the beam, including material properties, dimensions, mass density, damping ratio, and tip payload.
Mesh Spacing and Time Step Computation:

Calculates the spatial and temporal discretization parameters (dx and dt) for the finite difference method.
Initial Conditions:

Initializes the displacement w and disturbance d matrices and applies initial conditions reflecting the initial deflection due to the tip payload.
Frequency and Damping Coefficients Calculation:

Computes natural frequencies (w1, w2, ..., w6) and damping coefficient (ra) based on the physical parameters of the beam.
Boundary Disturbance and External Forces:

Defines the boundary disturbance d and external force distribution f along the beam.
Main Loop for Time-stepping:

Implements the finite difference method to update the displacement of the beam at each time step.
Applies boundary conditions at the beam's tip.
Visualization:

Plots the beam displacement over time using a 3D mesh plot.
Animates the displacement of the beam over time in a 2D plot.



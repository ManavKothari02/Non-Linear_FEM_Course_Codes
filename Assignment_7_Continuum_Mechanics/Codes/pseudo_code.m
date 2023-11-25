%% Pseudo Code
% Name: Manav Kothari
%UIN: 133008008

%% Assignment 7 

Inputs

Call MeshR function to generate mesh
Call PlotMesh function to plot the mesh

Define DOF_NOD 

Calculate the components of C

Initialize GCU, GPU, and GLS

for NL = 1:NLS % Load Step for loop
    Calculate the load for that particular step
    while iter<=ITERMAX && convergence == 0
        Initialize GLK and GLF

        for N = 1:NEM
            Calculate ELS and ELXY for every element

            Call ELEMATRICS2D function to calculate ELK and ELF

            Assemble the ELK and ELF matrix into GLK and GLF

        end

        Call CONTBCS function to impose boundary conditions

        Calculate iterative solution DELU
        % Note, the formulation of this method inherents NI algorithm. 
        % That means, this method already imposes NI

        for I = 1:NNM
            Update the solution vector and nodal coordinates
        end

        if MODEL ~= 1 
            Perform iterative method to solve for the solution 
        else
            Do not perform iteration and calculate solution based on 1st iteration itself
        end
    end % End of iterative loop

    %% Post processing of results

    if IGRAD ~= 0 %IGRAD = 0 means don't calculate stresses
        for I = 1:NPE
            Calculate ELXY and ELS from updated GLXY and GLS for the required element
            Call Stress2D function to get all the required stresses and strains.
        end
    end
end % End of load step loop

% --------------------------------------------------------------------------------------------- %

function [ELK,ELF] = ELEMATRICS2D(NDF,NPE,ELXY,ELS,NGPF,C,thick,F,LFORM)


Initialize ELK and ELF

for NI = 1:NGPF
    for NJ = 1:NGPF % Perform full integration
        Call INTERPLN2D function for SFL, GDSFL, and JAC

        Calculate Green Strain and 2nd Piola-Kirchhoff stress for Total lagrange
        OR
        Calculate Euler Strain and Cauchy stress for Updated Lagrange

        for I = 1:NPE
            Calculate components of ELF matrix
            for J = 1:NPE
                calculate components of ELK matrix
            end
        end
    end
end
end

% ---------------------------------------------------------------------------------------- %

function SS = STRESS2D(NDF,ELXY,ELS,LGP,NPE,C)

for NI = 1:LGP
    for NJ = 1:LGP

        Calculate XC, YC, U1X, U1Y, V1X, V1Y for the particular gauss point

        Calculate Euler strain and Cauchy Stress Tensor

        Update ELXY

        Calculate Green Strain and 2nd Piola Kirchhoff Stress tensor
    end
end

Export the required stress tensor and strains
end
% ---------------------------------------------------------------------------------------- %













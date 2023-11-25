%% Pseudo Code
% Do not run this

Give all Inputs

Call MeshR function to generate 2D Rectangular Mesh
Call PlotMesh function to get visual reprsentation of the generated mesh

Define DOF_NOD matrix 
% DOF_NOD matrix is but a matrix similar to NOD matrix, but rather than
% mapping global node number with element number, it maps global degree of
% freedom with element number

Initialize current and previous iteration solution (GCU & GPU)

%% Big Loop

for NL = 1:NLS % Initiate load step loop
    Define required constants like Rho and Mu
    while iter<=ITERMAX && convergence == 0 % Initiate iterative loop
        Initialize GLK and GLF matrices
        for N = 1:NEM
            Calculate ELU matrix from GPU & GCU & Acceleration parameter (GAMA1)
    
            Define ELXY from GLXY
    
            Call FLUIDMATRICS function to calculate ELK and ELF
    
            Perform Assembly of ELK & ELF into GLK & GLF
        end

        Call FLUIDBCS function to apply Essential and Natural BCs

        Calculate current iteration solution GCU

        Null out VSPV for NI after 1st iteration

        Calculate error and check for convergence

        iter = iter+1;
    end % END OF ITERATIVE LOOP 

    %% Post Processing of converged solution

    if IGRAD ~= 0 % Check if post-processing is required by the user or not
        for I = 1:NEM % Calculate for all elements
            Define ELXY and ELU

            Call STRESS2D function to calculate Pressure for all the gaussian points in the element

            Call Press_Calc script to print the required Pressure values
            % Note Press_Calc is question specific and needs to be changed
            % if we solve any other question (except 10.8.1 & 10.8.4)

        end
    end
end % END OF LOAD STEP LOOP

PRINT SOLUTIONS

% ---------------------------------------------------------------------- %

function [ELK,ELF] = FLUIDMATRICS(NDF,NPE,NONLIN,ELXY,ELU,RHOAMU,NGPF,GAMA2)

Initialize ELK and ELF

% Full integration 
for NI = 1:NGPF
    for NJ = 1:NGPF
        Calculate ELF and part of ELK matrix

        if NONLIN > 1
            Calculate TANG matrix as well
        end
    end
end

% Reduced Integration 
for NI = 1:NGPR
    for NJ = 1:NGPR
        Calcuate the penalty term in ELK
    end
end

if NONLIN > 1
    Calculate final TANG matrix and Residual matrix
end
end

% ------------------------------------------------------------------------%

function [XMAT,YMAT,SXMAT,SYMAT,SXYMAT,PRSMAT] = STRESS2D(ELXY,NPE,ELU,NGPR,GAMA2,MU)

% Perform Reduced integration
for NI = 1:NGPR
    for NJ = 1:NGPR
        Find the requried Gauss points where Pressure needs to be calculated in an element
        Calculate corresponding pressure.
    end
end







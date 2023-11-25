%% 2D 1DOF Finite Element Method code

clc, clear all, close all

%% Inputs
% Simulation Parameters
NONLIN = 2;
GAMA = 0;
ITERMAX = 10;
Epsilon = 0.001;
POSTP = 0; % POSTP == 1 if want to run postprocessing, else 0

% Differential Equation Coefficients
PDECOEFF.A11 = [0.2 0 0 0.04 0 0]; % A11 = [A10 A1x A1y A1u A1ux A1uy]
PDECOEFF.A22 = [0.2 0 0 0.04 0 0]; % A22 = [A20 A2x A2y A2u A2ux A2uy]
PDECOEFF.A00 = [0 0 0]; % A00 = [A00 A0x A0y]
PDECOEFF.F = [0 0 0]; % F = [F0 FX FY]

% Mesh Inputs
NX = 2; % Number of elements along x direction
NY = 2; % Number of elements along y direction
NEM = NX*NY; % Total number of elements in the mesh
NPE = 9;
if NPE<=4
    IEL = 1;
    NGPF = 2; % Number of Gauss Pts for Full integration
else
    IEL = 2;
    NGPF = 3; % Number of Gauss Pts for Full integration
end

% Rectangular Domain Dimesion
x_length = 0.2; % Length of rectangular domain in x axis
y_length = 0.1; % Length of rectangular domain in y axis

% Origin
X0 = 0; % X coord of element 1 node 1
Y0 = 0; % Y coord of element 1 node 1

% Essential Boundary Conditions
NSPV=13;
ISPV = [1 1; 6 1; 11 1; 16 1; 21 1; 5 1; 10 1; 15 1; 20 1; 25 1; 22 1; 23 1; 24 1];
VSPV=[500 500 500 500 500 300 300 300 300 300 45 0400 350];

% Natural Boundary Conditions
NSSV=3;
ISSV=[2 1; 3 1; 4 1];
VSSV=[0 0 0];

% Derived Constants
DX = (x_length/NX)*ones(NX,1); % x length of each element along x direction
DY = (y_length/NY)*ones(NY,1); % y length of each element along y direction

NDF = 1;
NNM = 25; % Total number of Nodes
NEQ = NNM * NDF; % Total number of equations in the problem
NN = NPE * NDF; % Total number of DOF per element

GCU = zeros(NEQ,1); % Global Solution from current iteration
GPU = zeros(NEQ,1); % Global Solution from previous iteration

%% Generating Mesh

[NOD,GLXY] = MeshR(IEL,NX,NY,NPE,DX,DY,X0,Y0);

%% Big Loop
iter = 1;
convergence = 0;

while iter<=ITERMAX && convergence == 0

    GLK = zeros(NEQ,NEQ); % Initialize global coefficient matrices and vectors
    GLF = zeros(NEQ,1); % Initialize global source vector

    % Loop to calculate Element Matrix for all Elements
    for N = 1:NEM
        if NONLIN > 0
            ELU = GAMA*GPU(NOD(N,:)) + (1-GAMA)*GCU(NOD(N,:));
        end
        ELXY = GLXY(NOD(N,:),:);
        
        [ELK,ELF] = ELEMATRC2D(NGPF,NPE,NONLIN,ELXY,ELU,PDECOEFF);
    
        GLK(NOD(N,:),NOD(N,:)) = GLK(NOD(N,:),NOD(N,:)) + ELK;
        GLF(NOD(N,:),1) = GLF(NOD(N,:),1) + ELF;
    end

    [GLK,GLF] = BNDRYUNSYM(NONLIN,NDF,NEQ,GLK,GLF,NSPV,NSSV,ISPV,ISSV,VSPV,VSSV);

    if(NONLIN==1)
        GPU = GCU;
        GCU = GLK\GLF;
    elseif(NONLIN==2)
        GPU = GCU;
        DELU    = GLK\GLF;
        GCU     = GCU + DELU;
    end
    
    % Null out VSPV for Newton's iteration after 1st iteration
    if iter == 1
        if NONLIN > 1
            VSPV(:) = 0;
        end
    end

    %GP1 = GCU;
    error(iter) = sqrt((sum((GPU - GCU).^2))/(sum(GCU.^2)));
    
    if (error(iter) <= Epsilon)
        convergence = 1;
        iter
        disp('Current Iteration solution ')
        disp(GCU)
    end
    iter
    disp('Current Iteration solution ')
    disp(GCU)
    iter = iter+1;

    % Post Processing of the solution

    if POSTP == 1
    disp('Flux Calc');
        for N = 1:NEM
            for I = 1:NPE
                NI = NOD(N,I);
                ELXY(I,1) = GLXY(NI,1);
                ELXY(I,2) = GLXY(NI,2);
                ELU(I) = GCU(NI);
            end
            NGPR = NGPF - 1;
            if convergence == 1
                POSTPROC(ELXY,NDF,NPE,ELU,NGPF,NGPR,PDECOEFF,iter,convergence);
                break;
            end
        end
    end


end
if convergence ~= 1
    fprintf('Error: Did Not Converge \n');
end








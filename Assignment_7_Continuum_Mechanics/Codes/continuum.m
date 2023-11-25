%% Nonlinear FEM course code 
%Assignment 7 - Continuum Mechanics
% Manav Kothari
% 133008008
clc, clear all, close all

%% Inputs

% Simulation Inputs
PLANE = 0; % PLANE  = 1 - Plane Strain ; else Plane Stress
ITERMAX = 25;
Epsilon = 0.001;
GAMA = 0.5;
LFORM = 2; % LFORM  == 1 - TOTAL LAGRANGE, ELSE UPDATED LAGRANGE
IGRAD = 1; % IGRAD = 1 - PERFORM POST-PROCESSING
MODEL = 1; %  MODEL = 1 - NO ITERATION, ELSE ITERATIONS
% DOMAIN Inputs
X0 = 0;
Y0 = 0;

x_length = 10;
y_length = 1;
thick = 0.1;

% Mesh Inputs
NX = 5;
NY = 1;
NPE = 8;
NDF = 2;

NEM = NX*NY;

DX = (x_length/NX)*ones(NX,1);
DY = (y_length/NY)*ones(NY,1);

if NPE <= 4
    IEL = 1;
    NGPF = 2;
else
    IEL = 2;
    NGPF = 3;
end
NGPF = 2;
% NGPR = NGPF - 1;

% Material Properties
E = [1 1]*1.2e7;            % E = [E1 E2]
G = [1 1 1]*4.6154e6;       % G = [G12 G13 G23]
v = [0.3 0.3 0.3];          % v = [v12 v21 v23]

% Loading Condition
NLS = 32;
DPx = zeros(NLS,1);
DPy = 500*ones(NLS,1);
F = [0 0]; % F = [Fx Fy]

% Essential BCs
NSPV = 6;
ISPV = [1,2; 1,1; 12,2; 12 ,1; 18,2; 18,1];
VSPV = [0 0 0 0 0 0];
% NSPV = 6;
% ISPV = [1 1; 1 2; 12 1; 12 2; 23 1; 23 2];
% VSPV = [0 0 0 0 0 0];

% Natural BCs
NSSV = 20;
ISSV = [19 2; 20 2; 21 2; 22 2; 23 2; 24 2; 25 2; 26 2; 27 2; 28 2; 2 2; 3 2; 4 2; 5 2; 6 2; 7 2; 8 2; 9 2; 10 2; 11 2];
VSSV = [0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.016667 0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.016667];
% NSSV = 20;
% ISSV = [24 2; 25 2; 26 2; 27 2; 28 2; 29 2; 30 2; 31 2; 32 2; 33 2; 2 2; 3 2; 4 2; 5 2; 6 2; 7 2; 8 2; 9 2; 10 2; 11 2];
% VSSV = [0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.016667 0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.033333 0.066667 0.016667]';

%% Generating Mesh
[NOD, GLXY, NNM] = MeshR(IEL,NX,NY,NPE,DX,DY,X0,Y0);
PlotMesh(NNM,GLXY)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------- %

%% Defining DOF_NOD
DOF_NOD = zeros(NEM,NDF*NPE);

if NPE == 4 % For linear element
    for i = 1:NEM
        DOF_NOD(i,:) = [ (NDF*NOD(i,1)-NDF+1):NDF*NOD(i,1) (NDF*NOD(i,2)-NDF+1):NDF*NOD(i,2) (NDF*NOD(i,3)-NDF+1):NDF*NOD(i,3) (NDF*NOD(i,4)-NDF+1):NDF*NOD(i,4)];
    end
elseif NPE == 9 % For 9 noded element
    for i = 1:NEM
        DOF_NOD(i,:) = [ (NDF*NOD(i,1)-NDF+1):NDF*NOD(i,1) (NDF*NOD(i,2)-NDF+1):NDF*NOD(i,2) (NDF*NOD(i,3)-NDF+1):NDF*NOD(i,3) (NDF*NOD(i,4)-NDF+1):NDF*NOD(i,4) (NDF*NOD(i,5)-NDF+1):NDF*NOD(i,5) (NDF*NOD(i,6)-NDF+1):NDF*NOD(i,6) (NDF*NOD(i,7)-NDF+1):NDF*NOD(i,7) (NDF*NOD(i,8)-NDF+1):NDF*NOD(i,8) (NDF*NOD(i,9)-NDF+1):NDF*NOD(i,9) ];
    end
elseif NPE == 8
    for i = 1:NEM
        DOF_NOD(i,:) = [ (NDF*NOD(i,1)-NDF+1):NDF*NOD(i,1) (NDF*NOD(i,2)-NDF+1):NDF*NOD(i,2) (NDF*NOD(i,3)-NDF+1):NDF*NOD(i,3) (NDF*NOD(i,4)-NDF+1):NDF*NOD(i,4) (NDF*NOD(i,5)-NDF+1):NDF*NOD(i,5) (NDF*NOD(i,6)-NDF+1):NDF*NOD(i,6) (NDF*NOD(i,7)-NDF+1):NDF*NOD(i,7) (NDF*NOD(i,8)-NDF+1):NDF*NOD(i,8) ];
    end
end
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------- %

%% Main Code

% Calculating C matrix
C = zeros(6);
if PLANE == 1
    denom = 1 - v(3) - 2*v(1)*v(2);
    C(1,1) = (E(1)*(1 - v(3)))/denom;
    C(1,2) = v(2)*E(1)/denom;
    C(2,1) = C(1,2);
    C(2,2) = E(2)*(1 - v(1)*v(2))/((1 + v(3))*denom);
    C(4,4) = G(2);
    C(5,5) = G(3);
    C(6,6) = G(1);
else
    C(1,1) = E(1)/(1 - v(1)*v(2));
    C(2,2) = E(2)/(1 - v(1)*v(2));
    C(1,2) = v(1)*C(2,2);
    C(2,1) = C(1,2);
    C(4,4) = G(2);
    C(5,5) = G(3);
    C(6,6) = G(1);
end

NN = NPE*NDF;
NEQ = NNM*NDF;

% Initial Guess
GCU = zeros(NEQ,1); % Current Iteration Solution
GPU = zeros(NEQ,1); % Previous Iteration Solution
GLS = zeros(NEQ,1); % Total Displacement 

MATRIX = zeros(NLS,5);

for NL = 1:NLS
    iter = 1;
    convergence = 0;
    F(1) = F(1) + DPx(NL);
    F(2) = F(2) + DPy(NL);
    
    while iter<=ITERMAX && convergence == 0
        GLK = zeros(NEQ);
        GLF = zeros(NEQ,1);

        for N = 1:NEM
            for I = 1:NPE
                NI = NOD(N,I);
                L = NI*NDF - 1;
                K = I*NDF - 1;
                ELS(K) = GLS(L);
                ELS(K+1) = GLS(L+1);
            end
            ELXY = GLXY(NOD(N,:),:);
            
            % Calculating Element Stiffness Matrix and Element force vector
            %[ELK,ELF] = CONTELE(NDF,NPE,ELXY,ELS,NGPF,C,thick,F,LFORM);
            [ELK,ELF] = ELEMATRICS2D(NDF,NPE,ELXY,ELS,NGPF,C,thick,F,LFORM);
            
            % Assembly
            GLK(DOF_NOD(N,:),DOF_NOD(N,:)) = GLK(DOF_NOD(N,:),DOF_NOD(N,:)) + ELK;
            GLF(DOF_NOD(N,:),1) = GLF(DOF_NOD(N,:),1) + ELF;
        end

        % Imposing BCs
        VSSVNL = VSSV*F(2); % !!!! QUESTION SPECIFIC !!!!
        [GLK,GLF] = CONTBCS(NDF,NEQ,GLK,GLF,NSPV,ISPV,VSPV,NSSV,ISSV,VSSVNL);

        % For this code, we will only be working with Newton's iteration.
        % This is because the derivation of stiffness matrix itself
        % involvesthe newton's iteration algorithm.

        % Calculating iteration solution
        DELU = GLK\GLF;
        
        GLF = DELU;
        % Update the solution vectors and nodal coordinates
        for I = 1:NNM
            L = (I-1)*NDF + 1;

            % Update the incremental displacement
            GPU(L) = GPU(L) + GLF(L);
            GPU(L+1) = GPU(L+1) + GLF(L+1);
            
            % Update total displacement
            GLS(L) = GLS(L) + GLF(L);
            GLS(L+1) = GLS(L+1) + GLF(L+1);
            
            % Update geometry in Updated Lagrangian Formulation
            if LFORM > 1
                GLXY(I,1) = GLXY(I,1) + GLF(L);
                GLXY(I,2) = GLXY(I,2) + GLF(L+1);
            end
        end
        if MODEL ~= 1
            SNORM = 0;      ENORM = 0;
            for I = 1:NEQ
                SNORM = SNORM + GLS(I)^2;
                ENORM = ENORM + GLF(I)^2;
            end
    
            TOL = sqrt(ENORM/SNORM);
            
            if TOL > Epsilon
                if iter > ITERMAX
                    break;
                else
                    iter = iter + 1;
                end
            else
                convergence = 1;
                break;
            end
        elseif MODEL == 1
            convergence = 1;
            break;
        end
    end % End of Iterative Loop

    if convergence ~= 1
        disp('Did not converge :(');
        break;
    end
    
    %disp("Load Step: " + NL);
    if NPE == 9
        MATRIX(NL,1) = thick*F(2);
        MATRIX(NL,4) = -GLS(43);
        MATRIX(NL,5) = GLS(44);
        MATRIX(NL,2) = GLXY(22,1); % X-COORD
        MATRIX(NL,3) = GLXY(22,2)-1; % Y-COORD
    elseif NPE == 8
        MATRIX(NL,1) = thick*F(2);
        MATRIX(NL,4) = -GLS(33);
        MATRIX(NL,5) = GLS(34);
        MATRIX(NL,2) = GLXY(17,1); % X-COORD
        MATRIX(NL,3) = GLXY(17,2)-1; % Y-COORD        
    end

    %% POST PROCESSING OF RESULTS
    ENUM = 1; % Element number
    if IGRAD ~= 0
         for I=1:NPE
            NI = NOD(ENUM,I);
            ELXY(I,1) = GLXY(NI,1);
            ELXY(I,2) = GLXY(NI,2);
            LI = (NI-1)*NDF;
            L  = (I-1)*NDF;
            for J=1:NDF
                LI = LI+1;
                L  = L+1;
                ELS(L) = GLS(LI);
            end
         end
        LGP = 2;
        SS = STRESS2D(NDF,ELXY,ELS,LGP,NPE,C);
        %SS.LOADSTEP = NL;
        STRESS_STRAIN(:,:,NL) = SS;

    end
end % End of Load Step loop

MATRIX2 = zeros(18,7);

for NL = 1:NLS
    MATRIX2(NL,1) = NL*50;
    MATRIX2(NL,2:7) = STRESS_STRAIN(1,:,NL)*10^5;
end












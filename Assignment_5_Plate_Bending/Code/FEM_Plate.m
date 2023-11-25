%% Assignment 5 - Plate Bending
% Manav Kothari
% UIN: 133008008

clc, clear all, clear all

%% Inputs
model = 2; % Not used for this code
% Mesh Inputs
NX = 4;
NY = 4;
NPE = 9;
NDF = 5;

NEM = NX*NY;

if NPE <= 4
    IEL = 1;
    NGPF = 2;
else
    IEL = 2;
    NGPF = 3;
end

% Domain Dimensions (Rectangular)
x_length = 5;
y_length = 5;
h = 1;

X0 = 0; % X coord of element 1 node 1
Y0 = 0; % Y coord of element 1 node 1

DX = (x_length/NX)*ones(NX,1); % x length of each element along x direction
DY = (y_length/NY)*ones(NY,1); % y length of each element along y direction

% Physical Parameters
E = [7.8e6,2.6e6]; % E = [E1 E2]
v12 = 0.25;
v21 = v12*E(2)/E(1);
v = [v12 v21];
G = 1.3*10^6*ones(1,3); % G = [G12 G13 G23]
Ks = 5/6;

% Loading conditions
% DelP = 0.05;
% P = 2;
Q0 = 32*5*((2.6e6*h^4)/(2*x_length)^4);
DelQ = 5*((2.6e6*h^4)/(2*x_length)^4);
DP=ones(1,Q0/DelQ);
NLS=max(size(DP));

NGPR = NGPF - 1;
NGPS = NGPF - 1;

% Simulation Parameters
NONLIN = 2; % NONLIN = 1 - DI ; NONLIN = 2 - NI
ITERMAX = 25;
Epsilon = 0.001;
GAMA = 0;
IGRAD = 1; % IGRAD = 1 - Calc stress


% Essential BCs
SS = 4; % SS = 1 - SS1 ; SS = 3 - SS3 ; SS = 4 - Clamped
ISPV_SS1 = [1 1; 1 2; 1 4; 1 5; 2 2; 2 5; 3 2; 3 5; 4 2; 4 5; 5 2; 5 5; 6 2; 6 5; 7 2; 7 5; 8 2; 8 5; 9 2; 9 3; 9 5; 10 1; 10 4; 19 1; 19 4; 28 1; 28 4; 37 1; 37 4; 46 1; 46 4; 55 1; 55 4; 64 1; 64 4; 73 1; 73 3; 73 4; 74 1; 74 3; 74 4; 75 1; 75 3; 75 4; 76 1; 76 3; 76 4; 77 1; 77 3; 77 4; 78 1; 78 3; 78 4; 79 1; 79 3; 79 4; 80 1; 80 3; 80 4; 81 1; 81 2; 81 3; 81 4; 81 5; 18 2; 18 3; 18 5; 27 2; 27 3; 27 5; 36 2; 36 3; 36 5; 45 2; 45 3; 45 5; 54 2; 54 3; 54 5; 63 2; 63 3; 63 5; 72 2; 72 3; 72 5];
load('ISPV_Clamped.mat');
load('ISPV_SS3.mat');

if SS == 1
    ISPV = ISPV_SS1;
elseif SS == 3
    ISPV = ispv_ss3';
elseif SS == 4
    ISPV = ispv_clamped';
end
NSPV = max(size(ISPV));

VSPV = zeros(NSPV,1);
% Natural BCs
NSSV = 0;
ISSV = [];
VSSV = [];


% -------------------------------------------------------------------------------------------- % 

%% Generate Mesh
[NOD,GLXY,NNM] = MeshR(IEL,NX,NY,NPE,DX,DY,X0,Y0);
PlotMesh(NNM,GLXY) % Plotting mesh

% -------------------------------------------------------------------------------------------- % 
%% Defining DOF_NOD
DOF_NOD = zeros(NEM,NPE*NDF);
if NPE == 4 % For linear element
    for i = 1:NEM
        DOF_NOD(i,:) = [ (NDF*NOD(i,1)-4):NDF*NOD(i,1) (NDF*NOD(i,2)-4):NDF*NOD(i,2) (NDF*NOD(i,3)-4):NDF*NOD(i,3) (NDF*NOD(i,4)-4):NDF*NOD(i,4)];
    end
elseif NPE == 9 % For 9 noded element
    for i = 1:NEM
        DOF_NOD(i,:) = [ (NDF*NOD(i,1)-4):NDF*NOD(i,1) (NDF*NOD(i,2)-4):NDF*NOD(i,2) (NDF*NOD(i,3)-4):NDF*NOD(i,3) (NDF*NOD(i,4)-4):NDF*NOD(i,4) (NDF*NOD(i,5)-4):NDF*NOD(i,5) (NDF*NOD(i,6)-4):NDF*NOD(i,6) (NDF*NOD(i,7)-4):NDF*NOD(i,7) (NDF*NOD(i,8)-4):NDF*NOD(i,8) (NDF*NOD(i,9)-4):NDF*NOD(i,9) ];
    end
elseif NPE == 8
    for i = 1:NEM
        DOF_NOD(i,:) = [ (NDF*NOD(i,1)-4):NDF*NOD(i,1) (NDF*NOD(i,2)-4):NDF*NOD(i,2) (NDF*NOD(i,3)-4):NDF*NOD(i,3) (NDF*NOD(i,4)-4):NDF*NOD(i,4) (NDF*NOD(i,5)-4):NDF*NOD(i,5) (NDF*NOD(i,6)-4):NDF*NOD(i,6) (NDF*NOD(i,7)-4):NDF*NOD(i,7) (NDF*NOD(i,8)-4):NDF*NOD(i,8) ];
    end
end
% -------------------------------------------------------------------------------------------- % 
%% Main Loop
NEQ = NNM*NDF;

GPU = zeros(NEQ,1); % Previous iteration solution
GCU = zeros(NEQ,1); % Current iteration solution

f0 = 0;
Center_Def = zeros(NLS,1);

for NL = 1:NLS
    iter = 1;
    f0 = f0 + DP(NL);
    q0 = f0*DelQ;
    convergence = 0;

    while iter<=ITERMAX && convergence == 0
        GLK = zeros(NEQ);
        GLF = zeros(NEQ,1);

        % Loop to calculate element matrix and assembling the matrix
        for N = 1:NEM
            if NONLIN > 0
                ELU = GAMA*GPU(DOF_NOD(N,:)) + (1-GAMA)*GCU(DOF_NOD(N,:));
            end

            ELXY = GLXY(NOD(N,:),:);
            
            [ELK,ELF] = FSDTPLTS(NGPF,NGPR,NGPS,NPE,NDF,NONLIN,ELXY,ELU,q0,E,v,G,Ks,h);

            GLK(DOF_NOD(N,:),DOF_NOD(N,:)) = GLK(DOF_NOD(N,:),DOF_NOD(N,:)) + ELK;
            GLF(DOF_NOD(N,:),1) = GLF(DOF_NOD(N,:),1) + ELF;
            
        end

        % Apply Boundary conditions
        [GLK,GLF] = BNDRYUNSYM(NONLIN,NDF,NEQ,GLK,GLF,NSPV,NSSV,ISPV,ISSV,VSPV,VSSV);

        if NONLIN == 1
            GPU = GCU;
            GCU = GLK\GLF;
        elseif NONLIN == 2
            GPU = GCU;
            DELU = GLK\GLF;
            GCU = GCU + DELU;
        end

        % Null out VSPV for NI after 1st iteration
        if iter == 1
            if NONLIN > 1
                VSPV(:) = 0;
            end
        end

        error(iter) = sqrt((sum((GPU - GCU).^2))/(sum(GCU.^2)));
        
        if error(iter) <= Epsilon
            convergence = 1;
        end
        iter = iter+1;
    end
    Center_Def(NL,1) = GCU(3);
    if convergence ~= 1
        disp('Error: Did not converge \n');
        break;
    end

    if convergence == 1    
        if IGRAD ~= 0
            disp('Stress:');
            for N = 1:NEM
                for I = 1:NPE
                    NI = NOD(N,I);
                    ELXY(I,1) = GLXY(NI,1);
                    ELXY(I,2) = GLXY(NI,2);
                    LI = (NI-1)*NDF;
                    L = (I-1)*NDF;
                    for J = 1:NDF
                        LI = LI + 1;
                        L = L + 1;
                        ELU(L) = GCU(LI);
                    end
                end
                ARRAY = [E(1) E(2) v12 h G(1) G(2) G(3)];
                Ele = 9;
                STRESS = POSTPROC2D(ELXY,NDF,NPE,ELU,NGPR,ARRAY);
                if Ele == N
                    disp(STRESS);
                    Stress(NL,1) = STRESS(1);
                    Stress(NL,3) = STRESS(3);
                    Stress(NL,2) = STRESS(2);
                    Stress(NL,5) = STRESS(4);
                    Stress(NL,4) = STRESS(5);
                end
            end
        end
    end
end







        










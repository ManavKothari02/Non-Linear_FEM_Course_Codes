%% Assignment 4 Time dependent 2D DE
clc, clear all, close all
% Manav Kothari
% 133008008
% -------------------------------------------------%
%% Inputs

% Mesh Input
NX = 8;
NY = 8;
NPE = 4;
NDF = 1;

NEM = NX*NY;

if NPE <= 4
    IEL = 1;
    NGPF = 2;
else
    IEL = 2;
    NGPF = 3;
end

% Domain Dimensions (Rectangular)
x_length = 1;
y_length = 1;

X0 = 0; % X coord of element 1 node 1
Y0 = 0; % Y coord of element 1 node 1

DX = (x_length/NX)*ones(NX,1); % x length of each element along x direction
DY = (y_length/NY)*ones(NY,1); % y length of each element along y direction

% Time Simulation parameter
DT = 0.05;
NTime = 25; % Total number of Time steps

End_Time = DT*NTime;

% Differential Equation Parameters
ITEM = 1; % ITEM = 1 - Parabolic Eqn ; ITEM = 2 - Hyperbolic Eqn; ITEM = 0 - Static Problem

PDECOEFF.A11 = [1 0 0 0 0 0]; % Axx = [A0 A1x A1y A1u A1ux A2uy]
PDECOEFF.A22 = [1 0 0 0 0 0];
PDECOEFF.A00 = [0 0 0 0 0 0];
PDECOEFF.C = [1 0 0];        % C0 = [C0 Cx Cy]
PDECOEFF.F = [0 0 0];    % F = [F0 Fx Fy]

%DL = (F/NLS)*ones(NLS,1); % Assuming uniform increment of load

% Simulation Parameters
NONLIN = 1; % NONLIN = 1 - DI ; 2 - NI
ITMAX = 5;
Epsilon = 0.001;
NLS = 5; % Number of Load steps
alfa = 1; % From the Alfa family (time approximation)
GAMA = 0.5; % Time approximation of Hyperbolic equation

% Essential Boundary Conditions
NSPV = 17;
ISPV = [9 1; 18 1; 27 1; 36 1; 45 1; 54 1; 63 1; 72 1;81 1; 80 1; 79 1; 78 1; 77 1; 76 1; 75 1; 74 1; 73 1];
VSPV = zeros(NSPV,1);

% Natural Boundary Condition
NSSV = 15;
ISSV = [1 1; 2 1; 3 1; 4 1; 5 1; 6 1; 7 1; 8 1; 10 1; 19 1; 28 1; 37 1; 46 1; 55 1; 64 1];
VSSV = zeros(NSSV,1);

% Additional Input
%Check_Time = [0.05:0.05:1 1.25]; % Time steps at which you want your answer
Check_Time = 0; % If checking for all time step
if length(Check_Time) == 1
    Final_Value = zeros(NTime,2);
else
    Final_Value = zeros(length(Check_Time),2);
end

% -------------------------------------------------%

%% Generate Mesh
[NOD,GLXY,NNM] = MeshR(IEL,NX,NY,NPE,DX,DY,X0,Y0);
%PlotMesh(NNM,GLXY)

% -------------------------------------------------%

%% Big Loop

NEQ = NNM * NDF;
% Initializing solution vectors
GPU = zeros(NEQ,1); % Current Time Step and Current Iteration solution
GLU = zeros(NEQ,1); % Current Time Step and Previous Iteration solution
GLP = zeros(NEQ,1); % Previous Time Step solution
GLV = zeros(NEQ,1); % Velocity matrix
GLA = zeros(NEQ,1); % Accn matrix

f = 0;
k = 0;
if ITEM == 1
    if NSSV > 0
        VSSV = VSSV*DT;
    end
end

A1 = alfa*DT;
A2 = (1-alfa)*DT;
if GAMA ~= 0
    A3 = 2/(GAMA*DT^2);
    A4 = A3*DT;
    A5 = (1-GAMA)/GAMA;
    A6 = (2*alfa)/(GAMA*DT);
    A7 = (2*alfa)/GAMA;
    A8 = ((alfa/GAMA)-1)*DT;
end

A_Coeff = [A1 A2 A3 A4 A5 A6 A7 A8];

for NT = 1:NTime
    Time = NT*DT;
    iter = 0;
    convergence = 0;
    
    while (iter <= ITMAX) && (convergence == 0)
        GLK = zeros(NEQ);
        GLF = zeros(NEQ,1);
        
        for N = 1:NEM
            if (NONLIN > 0) || (ITEM > 0) 
                ELU = GLU(NOD(N,:)); % Transfer current solution
                if ITEM>0
                    ELU0 = GLP(NOD(N,:)); % Transfer of previous time step
                end
            end

            if ITEM == 2
                ELU1 = GLV(NOD(N,:)); % Transfer of previous first time derivative
                ELU2 = GLA(NOD(N,:)); % Transfer of previous second time derivative
            else
                ELU1 = zeros(length(ELU),1);
                ELU2 = ELU1;
            end

            ELXY = GLXY(NOD(N,:),:);
            
            [ELK,ELM,ELF] = ELEMATRCS2D_Time(NGPF,NPE,NDF,NONLIN,ELXY,ELU,ELU0,ELU1,ELU2,PDECOEFF,ITEM,alfa,GAMA,DT,A_Coeff);
            
            % Assembly
            GLK(NOD(N,:),NOD(N,:)) = GLK(NOD(N,:),NOD(N,:)) + ELK;
            GLF(NOD(N,:),1) = GLF(NOD(N,:),1) + ELF;
        end

        [GLK,GLF] = BNDRYUNSYM_Time(NONLIN,NDF,NEQ,GLK,GLF,NSPV,NSSV,ISPV,ISSV,VSPV,VSSV);
        
        Soln = GLK \ GLF;

        if NONLIN > 0
            for I = 1:NEQ
                GPU(I) = GLU(I);
                if NONLIN == 1 % Direct Iteration
                    GLU(I) = Soln(I);
                else
                    GLU(I) = GLU(I) + Soln(I);
                end
            end
        end
        iter = iter + 1;
        
        if iter == 1
            if NONLIN == 2
                VSPV(:) = 0;
            end
        end

        error(iter) = sqrt((sum((GLU - GPU).^2))/(sum(GLU.^2)));
        
        if (error(iter) <= Epsilon)
            convergence = 1;
        end
    end

    Time = NT*DT;

    T = [0.25 0.5 1 1.25];
    for I = 1:length(T)
        if abs(Time-T(I))<Epsilon*0.1
            flag=1;
        end
    end

    if length(Check_Time) == 1
        disp("Time = " + Time)
        disp(GLU(1))
        variable(k+1,2) = GLU(1);
        variable(k+1,1) = Time;
        k = k+1;

    else
        for I = 1:length(Check_Time)
            
            if abs(Check_Time(I) - Time)<(Epsilon*0.1)
                disp("Time = " + Time)
                disp(GLU(1)) 

                Final_Value(I,1) = Time;
                Final_Value(I,2) = GLU(1);
                break;
            end
        end
    end
   
    if ITEM > 0
        DIFF = 0;
        Soln = 0;
        for I = 1:NEQ
            Soln = Soln + GLU(I)*GLU(I);
            DIFF = DIFF + (GLU(I) - GLP(I))^2;
            if ITEM == 2
                GLF(I) = A3*(GLU(I) - GLP(I)) - A4*GLV(I) - A5*GLA(I);
                GLV(I) = GLV(I) + A1*GLF(I) + A2*GLA(I);
                GLA(I) = GLF(I);
            end
        end
    end
    
    GLP = GLU;

    if convergence ~= 1
        fprintf('Error: Did Not Converge \n');
    end

end





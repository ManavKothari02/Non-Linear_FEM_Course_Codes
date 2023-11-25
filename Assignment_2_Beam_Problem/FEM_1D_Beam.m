%% 1D Problems Non-Linear FEM Code for EBT and TBT Beam problems

clc, clear all, close all

%% Input
% Solver input
model = 2; % mode = 1 - 1D Problem; mode = 2 - EBT ; mode = 3 - TBT
NONLIN = 2; %NONLIN = 1 - Direct Iteration ; NONLIN = 2 - Newton's Iteration
GAMA = 0.25;
ITERMAX = 25;
Epsilon = 0.001;

%Mesh Input
NEM = 8;
IEL = 1;
NDF = 3; %DOF per node

% Domain Input
X0 = 0;
AL = 100;
GLU = zeros((NEM*IEL+1)*NDF,1);

% Material Properties
E = 3e7; % Young's modulus
v = 0.3; % Poisson Ration
A = 1; % Cross sectional Area
I = 1/12; % Moment of Inertia
G = E/(2*(1+v)); % Shear Modulus
Ks = 10/12; % Shear Correction factor

% Loading Condition
F = [0 0 0]; % F = [FX0 FX1 FX2]
Q = [1 0 0]; % Q = [QX0 QX1 QX2]

DP = 0.5*ones(10/0.5);
NLS = length(DP);
NGP = 2;
LGP = NGP-1;

% Essential Boundary Conditions
NSPV = 4;
ISPV = [1,1;1,2;1,3;9,1];
VSPV = [0;0;0;0];

% Natural Boundary Condition
NSSV = 0;
ISSV = [0];
VSSV = [0];

% Mixed Boundary Condition
NSMB = 1;
BETA0 = 250;
BETAU = 0;
UREF = 0;
ISMB = [9,2];

%% Actual Code

NPE = IEL + 1;
NNM = NEM*IEL + 1; %Number of nodes in the mesh

Result_w_L2 = zeros(length(DP),2);
Result_w_L = zeros(length(DP),2);

NEQ = NNM * NDF; %Basically the dimension of stiffness matrix or the number of variable in U vector

%Defining GLX
DL = (AL/NEM)/IEL; %Distance between 2 nodes
GLX = zeros(NNM,1);
for i = 1:NNM
    GLX(i) = X0 + (DL*(i-1));
end

NOD = zeros(NEM,NPE);
for i = 1:NEM
    if i == 1
        NOD(i,:) = i:i+IEL;
    else
        NOD(i,:) = NOD(i-1,end):NOD(i-1,end)+IEL;
    end
end

DOF_NOD = zeros(NEM,NPE*NDF);

for i = 1:NEM
    DOF_NOD(i,:) = [(NDF*NOD(i,1)-2):(NDF*NOD(i,1)) (NDF*NOD(i,2)-2):(NDF*NOD(i,2))]';
end

GP2 = zeros(length(GLU),1);

f = 0;
q = 0;

%% Big Loop
for j = 1:NLS
    
    if length(DP) ~=1
        f = f + DP(j); %j-th step of increasing load
    else
        f = DP;
    end
    iter = 1;
    convergence = 0;

    while (iter<=ITERMAX) && (convergence == 0)
        GLK = zeros(NEQ);
        GLF = zeros(NEQ,1);
        
        GP1 = GLU;
        GP1 = GAMA*GP2 + (1-GAMA)*GP1;

        ELU = zeros(NPE*NDF,1);
        ELX = zeros(NPE,1);

        for N = 1:NEM
            ELU = GP1(DOF_NOD(N,:));
            ELX = GLX(NOD(N,:));
            
            syms x;
            
            [ELK,ELF] = ELMATRCS1D(model, IEL, NONLIN, ELX, ELU, NGP, LGP, E, A, I, G, Ks, f, F, Q);
            % ELK = Element Stiffness Matrix
            % ELF = Element source Matrix

            %Assembly
            GLK(DOF_NOD(N,:),DOF_NOD(N,:)) = GLK(DOF_NOD(N,:),DOF_NOD(N,:)) + ELK;
            GLF(DOF_NOD(N,:),1) = GLF(DOF_NOD(N,:),1) + ELF;
        end
        
        [GLK,GLF] = BNDRYUNS1D(NONLIN,NDF,GLK,GLF,GLU,NSPV,ISPV,VSPV,NSSV,ISSV,VSSV,NSMB,ISMB,BETA0,BETAU,UREF);

      if(NONLIN==1)
        GP2 = GLU;
        GLU = GLK\GLF;
      elseif(NONLIN==2)
        GP2 = GLU;
        DELU    = GLK\GLF;
        GLU     = GLU + DELU;
      end
%         Sol = GLK\GLF;
%         if NONLIN == 1
%             GP2 = GP1;
%             GP1 = Sol;
%         elseif NONLIN == 2
%             GP2 = GP1;
%             GP1 = GP1 + Sol;
%         end
        %Error Calc
        GP1 = GLU;
        error(iter) = sqrt((sum((GP2 - GP1).^2))/(sum(GP1.^2)));
    
        if (error(iter) <= Epsilon)
            convergence = 1;
            break;
        end
        
        iter = iter+1;

    end
    if convergence ~= 1
        fprintf('Error: Did Not Converge \n');
        break
    end
    Q_step = f*Q
    iter
    GP1
    
    if rem(NEM,2) == 0
        Result_w_L2(j,1) = f;
        Result_w_L2(j,2) = GP1((NEM/2+1)*NDF - 1);
        Result_w_L(j,1) = f;
        Result_w_L(j,2) = GP1(end-1);
    end

end
if convergence ==1 && rem(NEM,2)==0
    Result_w_L2
    Result_w_L
    figure(1)
    plot(Result_w_L2(:,1),Result_w_L2(:,2));
    xlabel('Load Steps');
    ylabel('Vertical deflection at L/2');
    title('Vertical Deflection at L/2 wrt Load Step');
    
    figure(2)
    plot(Result_w_L(:,1),Result_w_L(:,2));
    xlabel('Load Step');
    ylabel('Vertical Displacement at L');
    title('Vertical deflection at L vs Load Step')
end










%% 1D Problems Non-Linear FEM Code

clc, clear all, close all

%%Variable Definition 

% NNM: Number of nodes in the mesh
% NEM: Number of elements in the mesh
% NDF: Number of Degree of freedom
% AL: Domain Length
% NHBW: No of Half Bandwidth
% NBW: No of Bandwidth (2*NHBW)
% NPE: Node per element
% NOD: Connectivity Matrix
% IEL: Element Type (IEL = 1 - Linear ; IEL = 2 - Quadratic)
% GLK: Global coefficient matrix, [K]; xNote that the global source vector is stored in the last column of [GLK]
% GLX: Vector of Global coord of global nodes
% GLU: Vector of Global nodal values of Primary Variable
% GP1: r-1 iteration solution of PV
% GP2: r-2 iteration solution of PV
% NOD: Connectivity Matrix (related local node no to global node no)
% NONLIN: Flag for type of analysis
%         =0 for linear analysis
%         =1 for Direct Iteration method
%         >1 for Newton Iteration method 
% NLS: No of load steps
% ITMAX: Max iteration allowed
% EPS: Error Tolerance Specified
% GAMA: Acceleration parameter
% NEQ: No of equations in the mesh 
% NGP: Number of Gaussian Points
% ELX: Vector of global coord of the element nodes
% GJ: Jacobian of Transformation
% SFL: Shape interpolation funtion of Lagrange type
% GDSFL: Global drivative (wrt x) of SFL


%% Input Definition 

File = fopen('Input_File.inp');

for k = 1:29
    Current = fgetl(File);
    if (k >= 11)
        Current = split(Current);
        VALUE = str2num(char(Current(2)));
    else
        continue;
    end

    if k == 11
        NONLIN = VALUE;
    elseif k == 12
        X0 = VALUE;
    elseif k == 13
        AL = VALUE;
    elseif k == 14
        NEM = VALUE;
    elseif k == 15
        IEL = VALUE;
    elseif k == 16
        GAMA = VALUE;
    elseif k == 17
        ITERMAX = VALUE;
    elseif k == 18
        Epsilon = VALUE;
    elseif k == 19
        GP1 = VALUE;
    elseif k == 20
        A = VALUE;
    elseif k == 21
        B = VALUE;
    elseif k == 22
        C = VALUE;
    elseif k == 23
        F = VALUE;
    elseif k == 24
        NSPV = VALUE;
    elseif k == 25
        VSPV = VALUE;
    elseif k == 26
        ISPV = VALUE;   
    elseif k == 27
        NSSV = VALUE;
    elseif k == 28
        VSSV = VALUE;
    elseif k == 29
        ISSV = VALUE;
    end
end

%%Derived variables

NPE = IEL + 1;
NNM = NEM*IEL + 1; %Number of nodes in the mesh

NDF = 1; %DOF per node
NEQ = NNM * NDF; %Basically the dimension of stiffness matrix or the number of variable in U vector

NGP = IEL+1; %Number of Gaussian Quadrature points


%% Actual Code

%Defining GLX
DL = (AL/NEM)/IEL; %Distance between 2 nodes
GLX = zeros(NNM,1);
for i = 1:NNM
    GLX(i) = X0 + (DL*(i-1));
end

NOD = zeros(NEM,IEL+1);
for N = 1:NEM
    if N == 1
        NOD(N,:) = N:N+IEL;
    else
        NOD(N,:) = NOD(N-1,end):NOD(N-1,end)+IEL;
    end
end

% Big Loop
iter = 1;

convergence = 0;

GP1 = GP1'; % Making the guess vector a column vector

while (iter<=ITERMAX) && (convergence == 0)
    
    GLK = zeros(NEQ,NEQ); %Global Stiffness Matrix
    GLF = zeros(NEQ,1); %Global Force Matrix

    if iter == 1
        GP2 = zeros(NEQ,1);
    end

    GLU = zeros(NEQ,1);

    for i = 1:NEQ
        GLU(i) = GAMA*GP2(i) + (1-GAMA)*GP1(i);
    end

    ELU = zeros(NPE,1); %Primary Variable vector at element nodes
    ELX = zeros(NPE,1); %Global Coordinates of element nodes

    for N = 1:NEM
        for i = 1:NPE
            NI = NOD(N,i);
            ELU(i) = GLU(NI);
            ELX(i) = GLX(NI);
        end
        [ELK, ELF] = ELMATRCS1D(IEL, NONLIN, ELX, ELU, A, B, C, F);
        % ELK = Element Stiffness Matrix
        % ELF = Element source Matrix
        [GLK,GLF] = Assembly(ELK,ELF,NOD,N,GLK,GLF);
    end
    
    if NSPV>0
        for NP = 1:NSPV
            NB = (ISPV(NP,1)-1)*NDF + ISPV(NP,2);
            if NONLIN == 1 %For Direct Iteration Method
                GLF(NB) = VSPV(NP);
            elseif NONLIN == 2 %For Newton's Method
                GLF(NB) = 0;
            end
            GLK(NB,:) = 0;
            GLK(NB,NB) = 1;
        end

        if NSSV>0
            for i = 1:NSSV
                NB = (ISSV(i,1)-1)*NDF+ISSV(i,2);
                GLF(NB)=GLF(NB)+VSSV(i);
            end
        end
    end
    
    %Solution
    Sol = linsolve(GLK,GLF);
    
    if NONLIN == 1
        GP2 = GP1;
        GP1 = Sol;
    elseif NONLIN == 2
        GP2 = GP1;
        GP1 = GP1 + Sol;
    end
    
    %Error Calc
    error(iter) = sqrt((sum((GP2 - GP1).^2))/(sum(GP1.^2)));

    if (error(iter) <= Epsilon)
        convergence = 1;
    end
    iter
    iter = iter+1;
GP1

end

%% Analytical Solution
x_new = 0:0.01:1;
u = zeros(length(x_new),1);
for i = 1:length(u)
    u(i) = 1/(1+x_new(i));
end
ans = [GLX GP1]

figure(1)
plot([1:(iter-1)],error,'LineWidth',2);
xlabel('Iterations')
ylabel('Error');
title('Error vs iteration')

figure(2)
plot(x_new,u,GLX,GP1,'LineWidth',2);
xlabel('x');
ylabel('Primary Variable');
title('Primary variable vs x')
legend('Analytical','Simulation');  




%% -------------- Element Matrixed Function --------------------

% Elemental Stiffness and source matrix formation function
function [ELK,ELF] = ELMATRCS1D(IEL,NONLIN,ELX,ELU, A, B, C, F)
AX0 = A(1);     BX0 = B(1);     CX0 = C(1);     FX0 = F(1);
AX1 = A(2);     BX1 = B(2);     CX1 = C(2);     FX1 = F(2);
AU1 = A(3);     BU1 = B(3);     CU1 = C(3);     FX2 = F(3);
AU2 = A(4);     BU2 = B(4);     CU2 = C(4);     
AUX1 = A(5);    BUX1 = B(5);    CUX1 = C(5);
AUX2 = A(6);    BUX2 = B(6);    CUX2 = C(6);

NGP = IEL+1;

[GAUSPT, GAUSWT] = Gaus_int(NGP);

NPE = IEL + 1;
EL = ELX(end) - ELX(1); %Element Length

ELK = zeros(NPE);
ELF = zeros(NPE,1);

if NONLIN == 2
    T = zeros(NPE); %K Tangent
end

for NI = 1:NGP
    XI = GAUSPT(NI);
    [SFL, GDSFL,GJ] = INTERPLN1D(ELX,IEL,XI);
    X = ELX(1)+0.5*(1+XI)*EL;
    Const = GJ*GAUSWT(NI);

    U = 0;
    dU = 0;
    for i = 1:NPE
        U = U+SFL(i)*ELU(i);
        dU = dU+GDSFL(i)*ELU(i);
    end
    
    A = AX0 + AX1*X + AU1*U + AU2*U^2 + AUX1*dU + AUX2*dU^2;
    B = BX0 + BX1*X + BU1*U + BU2*U^2 + BUX1*dU + BUX2*dU^2;
    C = CX0 + CX1*X + CU1*U + CU2*U^2 + CUX1*dU + CUX2*dU^2;
    F = FX0 + FX1*X + FX2*X^2;

    if NONLIN == 2
        AXT1 = (AU1+2*AU2*U)*dU;
        AXT2 = (AUX1+2*AUX2*dU)*dU;
        BXT1 = (BU1+2*BU2*U)*dU;
        BXT2 = (BUX1+2*BUX2*dU)*dU;
        CXT1 = (CU1 + 2*CU2*U)*dU;
        CXT2 = (CUX1 + 2*CUX2*dU)*dU;
    end
    for i = 1:NPE
        ELF(i) = ELF(i) + F*SFL(i)*Const;
        for j = 1:NPE
            S00 = SFL(i)*SFL(j)*Const;
            S01 = SFL(i)*GDSFL(j)*Const;
            S10 = GDSFL(i)*SFL(j)*Const;
            S11 = GDSFL(i)*GDSFL(j)*Const;
            ELK(i,j) = ELK(i,j) + A*S11 + B*S01 + C*S00;
            if NONLIN == 2
                T(i,j) = T(i,j) + AXT1*S10 + AXT2*S11 + BXT1*S00 + BXT2*S01 + CXT1*S00 + CXT2*S01;
            end
        end
    end
end

if NONLIN == 2
    for i = 1:NPE
        for j = 1:NPE
            ELF(i) = ELF(i) - ELK(i,j)*ELU(j); % for Newton method, ELF is Residual and not soruce array
        end
    end
    ELK = ELK + T;
end
end


%% ------------ Shape function and its derivative -------------------------
function [SFL,GDSFL,GJ] = INTERPLN1D(ELX,IEL,XI)
syms x;

NPE = IEL+1;

if (IEL == 1) %ie linear element
    SFL = [0.5*(1-x) 0.5*(1+x)]; 
    DSFL = diff(SFL);
    x = XI;
    SFL = eval(SFL);
    DSFL = eval(DSFL);
elseif (IEL == 2) %ie quadratic element 
    SFL = [-0.5*x*(1-x) 1-x^2 0.5*x*(1+x)];
    DSFL = diff(SFL);
    x = XI;
    SFL = eval(SFL);
    DSFL = eval(DSFL);
end
GJ = 0; %Jacobian
for i = 1:NPE
    GJ = GJ + DSFL(i)*ELX(i);
end
GDSFL = DSFL/GJ;
end

%% -------------------- Gaussian points --------------------------------
function [x,w] = Gaus_int(n)

%Works for n = 2 to 5

if n == 2
    x = [-0.5773502691896257 0.5773502691896257]';
    w = [1 1]';
elseif n==3
    x = [0 -0.7745966692414834 0.7745966692414834]';
    w = [0.8888888888888888 0.5555555555555556 0.5555555555555556]';
elseif n==4
    x = [-0.3399810435848563 0.3399810435848563 -0.8611363115940526 0.8611363115940526]';
    w = [0.6521451548625461 0.6521451548625461 0.3478548451374538 0.3478548451374538]';
elseif n==5
    x = [0 -0.5384693101056831 0.5384693101056831 -0.9061798459386640 0.9061798459386640]';
    w = [0.5688888888888889 0.4786286704993665 0.4786286704993665 0.2369268850561891 0.2369268850561891]';
end
end

%% ------------------------ Assembly function -----------------------------------
function [GLK,GLF] = Assembly(ELK,ELF,NOD,N,GLK,GLF)
% N is the node number under consideration
% NEQ is total number of equations, size of GLK and GLF

GLK(NOD(N,:),NOD(N,:)) = GLK(NOD(N,:),NOD(N,:)) + ELK;
GLF(NOD(N,:),1) = GLF(NOD(N,:),1) + ELF;

end
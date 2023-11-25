%% FEM Viscous flow
% Manav Kothari
% UIN: 133008008

clc, clear all, close all

%% Inputs
Problem = 2;
% Simulation Inputs
NONLIN = 1;
ITERMAX = 25;
Epsilon = 0.001;
GAMA1 = 0.5; % Acceleration Parameter
GAMA2 = 10^8; % Penalty constant
IGRAD = 1; % IGRAD = 1 - perform post-processing ; IGRAD = 0 - No post processing

% Domain Inputs
X0 = 0;
Y0 = 0;
if Problem == 1
    x_length = 6;
    y_length = 2;
elseif Problem == 2
    x_length = 1;
    y_length = 1;
end

% Mesh Inputs
NX = 8; % Number of elements in x-direction
NY = 10; % Number of elements in y-direction
NPE = 9;

if Problem == 1
    if NX == 10 && NY == 6
        DX = [1 1 1 1 0.5 0.5 0.25 0.25 0.25 0.25]; % x length of each element along x direction
        DY = [0.25 0.25 0.5 0.5 0.25 0.25]; % y length of each element along y direction
    elseif NX == 5 && NY == 3
        DX = [2 2 1 0.5 0.5];
        DY = [0.5 1 0.5];
    end
elseif Problem == 2
    if NX == NY
        DX = (x_length/NX)*ones(1,8);
        DY = (x_length/NY)*ones(1,8);
    elseif NX == 8 && NY == 10
        DX = 0.125*ones(1,NX);
        DY = [0.125 0.125 0.125 0.125 0.125 0.125 0.0625 0.0625 0.0625 0.0625];
    elseif NX == 16 && NY == 20
        DX = (x_length/NX)*ones(1,NX);
        DY = [0.0625 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625 0.03125 0.03125 0.03125 0.03125 0.03125 0.03125 0.03125 0.03125];
    end
end

NDF = 2;
NEM = NX*NY;
if NPE <= 4
    IEL = 1;
    NGPF = 2;
else
    IEL = 2;
    NGPF = 3;
end
NGPR = NGPF - 1;

% Loading condition input
if Problem == 1
    DP = 1;
elseif Problem == 2
    if NY < 9
        DP = 250:250:750;
    else
        DP = 0:500:1000;
    end
end
F = [0 0]; % F = [FX FY]
MU = 1;

% Essenstial Boundary Conditions
if Problem == 1    
    ISPV = [1	77	76	75	74	73	72	71	70	69	68	67	56	45	34	23	12	1	2	3	4	5	6	7	8	9	10	11	77	76	75	74	73	72	71	70	69	68	67
            1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2]';
    VSPV = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1]';
    NSPV = length(VSPV);
elseif Problem == 2
    if NX == NY
        ISPV = [1	2	3	4	5	6	7	8	9	18	27	36	45	54	63	72	81	80	79	78	77	76	75	74	73	64	55	46	37	28	19	10	1	2	3	4	5	6	7	8	9	18	27	36	45	54	63	72	81	80	79	78	77	76	75	74	73	64	55	46	37	28	19	10
                1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2]';
        VSPV = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0]';
        NSPV = length(VSPV);
        Mat = zeros(7,4);
    elseif NX ~= NY 
        ISPV = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	34	51	68	85	102	119	136	153	170	187	204	221	238	255	272	289	306	323	340	357	356	355	354	353	352	351	350	349	348	347	346	345	344	343	342	341	324	307	290	273	256	239	222	205	188	171	154	137	120	103	86	69	52	35	18	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	34	51	68	85	102	119	136	153	170	187	204	221	238	255	272	289	306	323	340	357	356	355	354	353	352	351	350	349	348	347	346	345	344	343	342	341	324	307	290	273	256	239	222	205	188	171	154	137	120	103	86	69	52	35	18
                1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2]';
        VSPV = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0]';
        NSPV = length(VSPV);
        Mat = zeros(20,4);  
    end
end

% Natual BCs
if Problem == 1    
    ISSV = [1	2	3	4	5	6	7	8	9	10	11	22	33	44	55	66	77	1	11	22	33	44	55	66	77	67	56	45	34	23	12
            1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2]';
    VSSV = zeros(max(size(ISSV)),1);
    NSSV = length(VSSV);
elseif Problem == 2
    NSSV = 0;
    ISSV = [];
    VSSV = [];
end

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------- %

%% Generating Mesh
[NOD, GLXY, NNM] = MeshR(IEL,NX,NY,NPE,DX,DY,X0,Y0);
%PlotMesh(NNM,GLXY)
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

%% Main code

NN = NPE*NDF;
NEQ = NNM*NDF;

% Initial guess
GCU = zeros(NEQ,1);
GPU = zeros(NEQ,1);

NLS = max(size(DP));
F0 = 0;



for NL = 1:NLS
    iter = 1;
    convergence = 0;

    REYNOLDS = DP(NL);  VC = 1;     LC = 1;
    RHOAMU(1) = REYNOLDS*(VC*LC/MU);
    RHOAMU(2) = MU;

    while iter<=ITERMAX && convergence == 0
        GLK = zeros(NEQ);
        GLF = zeros(NEQ,1);

        for N = 1:NEM
            if NONLIN > 0
                ELU = GAMA1*GPU(DOF_NOD(N,:)) + (1-GAMA1)*GCU(DOF_NOD(N,:));
            elseif NONLIN == 0
                ELU = 0;
            end

            ELXY = GLXY(NOD(N,:),:);
            
            % Calculating Element stiffness and source matrix
            [ELK,ELF] = FLUIDMATRICS(NDF,NPE,NONLIN,ELXY,ELU,RHOAMU,NGPF,GAMA2);
            
            % Assembly
            GLK(DOF_NOD(N,:),DOF_NOD(N,:)) = GLK(DOF_NOD(N,:),DOF_NOD(N,:)) + ELK;
            GLF(DOF_NOD(N,:),1) = GLF(DOF_NOD(N,:),1) + ELF;
        end
        % Imposing BCs
        [GLK,GLF] = FLUIDBCS(NONLIN,NDF,NEQ,GLK,GLF,NSPV,ISPV,VSPV,NSSV,ISSV,VSSV);

        if NONLIN == 1
            GPU = GCU;
            GCU = GLK\GLF;
        elseif NONLIN == 2
            GPU = GCU;
            DELU = GLK\GLF;
            GCU = GCU + DELU;
        elseif NONLIN == 0 % Linear mode
            GCU = GLK\GLF;
            disp('Linear Mode');
            convergence = 1;
            iter = iter + 1;
            break;
        end
        
        % Null out VSPV for NI after 1st iteration
        if iter == 1
            if NONLIN > 1
                VSPV(:) = 0;
            end
        end

        error = sqrt((sum((GPU - GCU).^2))/(sum(GCU.^2)));
    
        if error <= Epsilon
            convergence = 1;
        end
        iter = iter+1;
    end % End of While iterative loop

    if convergence ~= 1
        disp('Error: Did not converge :(');
        break;
    end
    
    if Problem == 2
        if NY < 9
            A = 14:9:77;
            for I = 1:7
                Mat(I,1) = GLXY(A(I),2);
                Mat(I,NL+1) = GCU(2*A(I)-1);
            end
        else
            A = 26:17:349;
            for I = 1:20
                Mat(I,1) = GLXY(A(I),2);
                Mat(I,NL+1) = GCU(2*A(I)-1);
            end
        end
    end

    %% Post processing of the solution
    if IGRAD ~= 0
        for N = 1:NEM
            ELXY = GLXY(NOD(N,:),:);
            ELU = GCU(DOF_NOD(N,:));

            [X,Y,SX,SY,SXY,PRS] = STRESS2D(ELXY,NPE,ELU,NGPR,GAMA2,MU);
            
            MATRIX{N,1,NL} = X;
            MATRIX{N,2,NL} = Y;
            MATRIX{N,3,NL} = SX;
            MATRIX{N,4,NL} = SY;
            MATRIX{N,5,NL} = SXY;
            MATRIX{N,6,NL} = PRS;
            
            Pressure_Calc
        end
    end

end % End of Load step loop

%% Print Solution
fprintf('iter: %d \n',iter);
fprintf('Load Step: %d \n', NL);

fprintf('Solution: \n');
Disp_V = [];

if Problem == 1
    Mat1 = zeros(7,2);
    Mat2 = zeros(7,2);

    II = 3;
    for I = 1:10
        Disp_V(length(Disp_V)+1) = GCU(II);
        II = II+2;
    end
    Disp_V = Disp_V'
    % Vx at x = 4
    A = 5:11:71;
    for I = 1:length(A)
        Mat1(I,2) = GCU(2*A(I)-1);
        Mat1(I,1) = GLXY(A(I),2);
    end
    disp('Vx at x = 4');
    Mat1

    % Vx at x = 6
    A = 11:11:77;
    for I = 1:length(A)
        Mat2(I,2) = GCU(2*A(I)-1);
        Mat2(I,1) = GLXY(A(I),2);
    end
    disp('Vx at x = 6');
    Mat2

elseif Problem == 2
    disp('vx(0.5,y)');
    Mat
end









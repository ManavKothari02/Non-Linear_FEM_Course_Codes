function [NOD,GLXY] = MeshR(IEL,NX,NY,NPE,DX,DY,X0,Y0)

% NOD = zeros(MAXNEM,9); % Nodal Connectivity
% GLXY = zeros(MAXNNM,2); % Global Coords of nodes
% DX = zeros(MAXNX,1); % x coords of nodes
% DY = zeros(MAXNY,1); % y coords of nodes

NEX1 = NX+1;
NEY1 = NY+1;
NXX = IEL*NX;
NYY = IEL*NY;
NXX1 = NXX+1; % Total number of nodes along x direction
NYY1 = NYY+1; % Total number of nodes along y direction

NEM = NX*NY;
NNM = NXX1*NYY1; % Total number of nodes 

if NPE == 8
    NNM = NNM - NEM;
end

% Rectangular Element
K0 = 0;
if NPE == 9
    K0 = 1;
end

NOD(1,1) = 1;
NOD(1,2) = IEL+1;
NOD(1,3) = NXX1 + (IEL-1)*NEX1 + IEL + 1;

if NPE == 9
    NOD(1,3) = 4*NX + 5;
end

NOD(1,4) = NOD(1,3) - IEL;

if NPE>4
    NOD(1,5) = 2;
    NOD(1,6) = NXX1 + (NPE-6);
    NOD(1,7) = NOD(1,3) - 1;
    NOD(1,8) = NXX1 + 1;

    if NPE == 9
        NOD(1,9) = NXX1+2;
    end
end

if NY>1
    M = 1;
    for N = 2:NY
        L = (N-1)*NX + 1;
        for I = 1:NPE
            NOD(L,I) = NOD(M,I) + NXX1 + (IEL-1)*NEX1 + K0*NX;
        end
        M = L;
    end
end

if NX>1
    for NI = 2:NX
        for I = 1:NPE
            K1 = IEL;
            if (I == 6) || (I == 8)
                K1 = 1+K0;
            end
            NOD(NI,I) = NOD(NI-1,I) + K1;
        end
        M = NI;

        for NJ = 2:NY
            L = (NJ-1)*NX + NI;
            for J = 1:NPE
                NOD(L,J) = NOD(M,J) + NXX1 + (IEL-1)*NEX1 + K0*NX;
            end
            M = L;
        end
    end
end

%% Generate Global Coordinates of the nodes

DX(NEX1) = 0;
DY(NEY1) = 0;

XC = X0;
YC = Y0;

if NPE == 8
    for NI = 1:NEY1
        I = (NXX1 + NEX1)*(NI - 1) + 1;
        J = 2*NI - 1;
        GLXY(I,1) = XC;
        GLXY(I,2) = YC;
        for NJ = 1:NX
            DELX = 0.5*DX(NJ);
            I = I + 1;
            GLXY(I,1) = GLXY(I-1,1) + DELX;
            GLXY(I,2) = YC;
            I = I + 1;
            GLXY(I,1) = GLXY(I-1,1) + DELX;
            GLXY(I,2) = YC;
        end
        if NI<=NY
            I = I + 1;
            YC = YC + 0.5*DY(NI);
            GLXY(I,1) = XC;
            GLXY(I,2) = YC;
            for II = 1:NX
                I = I + 1;
                GLXY(I,1) = GLXY(I-1,1) + DX(II);
                GLXY(I,2) = YC;
            end
        end
    YC = YC + 0.5*DY(NI);
    end
else
    YC = Y0;
    for NI = 1:NEY1
        XC = X0;
        I = NXX1*IEL*(NI-1);
        for NJ = 1:NEX1
            I = I + 1;
            GLXY(I,1) = XC;
            GLXY(I,2) = YC;
            if NJ<NEX1
                if IEL == 2
                    I = I +1;
                    XC = XC + 0.5*DX(NJ);
                    GLXY(I,1) = XC;
                    GLXY(I,2) = YC;
                end
            end
            XC = XC + DX(NJ)/IEL;
        end
        XC = X0;
        if IEL == 2
            YC = YC + 0.5*DY(NI);
            for NJ = 1:NEX1
                I = I + 1;
                GLXY(I,1) = XC;
                GLXY(I,2) = YC;
                if NJ<NEX1
                    I = I + 1;
                    XC = XC + 0.5*DX(NJ);
                    GLXY(I,1) = XC;
                    GLXY(I,2) = YC;
                end
                XC = XC + 0.5*DX(NJ);
            end
        end
        YC = YC + DY(NI)/IEL;
    end
end

end
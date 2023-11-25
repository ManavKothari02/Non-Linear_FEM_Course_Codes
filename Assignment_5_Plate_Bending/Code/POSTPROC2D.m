function STRESS = POSTPROC2D(ELXY,NDF,NPE,ELU,NGPR,ARRAY)

[GAUSPT,GAUSWT] = Gaus_int(NGPR);

for NI = 1:NGPR
    for NJ = 1:NGPR
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        XC = 0;     YC = 0;

        for I = 1:NPE
            XC = XC + SFL(I)*ELXY(I,1);
            YC = YC + SFL(I)*ELXY(I,2);
        end

        
    end
end

E1 = ARRAY(1);      E2 = ARRAY(2);      v12 = ARRAY(3);
h = ARRAY(4);       G12 = ARRAY(5);     ZC = h/2;
G13 = ARRAY(6);     G23 = ARRAY(7);

SX = 0;     SY = 0;     DUX = 0;
DUY = 0;    DVX = 0;    DVY = 0;
DWX = 0;    DWY = 0;    DSXY = 0;
DSYX = 0;   DSXX = 0;   DSYY = 0;

for I = 1:NPE
    J = NDF*(I-1) + 1;
    J1 = J+1;
    J2 = J+2;
    K = J2+1;
    L = K+1;

    DUX  = DUX+GDSFL(1,I)*ELU(J);
    DUY  = DUY+GDSFL(2,I)*ELU(J);
    DVX  = DVX+GDSFL(1,I)*ELU(J1);
    DVY  = DVY+GDSFL(2,I)*ELU(J1);
    DWX  = DWX+GDSFL(1,I)*ELU(J2);
    DWY  = DWY+GDSFL(2,I)*ELU(J2);
    SX  = SX+SFL(I)*ELU(K);
    SY  = SY+SFL(I)*ELU(L);
    DSXX = DSXX+GDSFL(1,I)*ELU(K);
    DSXY = DSXY+GDSFL(2,I)*ELU(K);
    DSYX = DSYX+GDSFL(1,I)*ELU(L);
    DSYY = DSYY+GDSFL(2,I)*ELU(L);
end

v21 = v12*E2/E1;
EXX = DUX + 0.5*DWX^2+ZC*DSXX;
EYY = DVY + 0.5*DWY^2+ZC*DSYY;
GXY = DUY + DVX + DWX*DWY + ZC*(DSXY+DSYX);
Q11 = E1/(1-v12*v21);
Q12 = v12*E2/(1-v12*v21);
Q22 = E2/(1-v12*v21);
Q66 = G12;
Q = [Q11,Q12,0;Q12,Q22,0;0,0,Q66];
SIGMA = Q*[EXX;EYY;GXY];
SCF=5.0D0/6.0D0;
C44=SCF*G23*h;
C55=SCF*G13*h;
SGMXZ = 1.2*C55*(DWX+SX)/h;
SGMYZ = 1.2*C44*(DWY+SY)/h;
SIGMA = [SIGMA;SGMXZ;SGMYZ];
STRESS = SIGMA;


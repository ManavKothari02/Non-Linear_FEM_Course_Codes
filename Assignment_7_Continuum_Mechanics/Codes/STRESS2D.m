function SS = STRESS2D(NDF,ELXY,ELS,LGP,NPE,C)

[GAUSPT,GAUSWT] = Gaus_int(LGP);
COUNTER = 1;
CAUCHY_STRESS = zeros(LGP*LGP,3);
EULER_STRAIN = CAUCHY_STRESS;
PIOLA_STRESS = EULER_STRAIN;
GREEN_STRAIN = PIOLA_STRESS;

XY0_ARRAY = zeros(LGP*LGP,2);
XYC_ARRAY = zeros(LGP*LGP,2);

for NI = 1:LGP
    for NJ = 1:LGP
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);

        U1X = 0;        U1Y = 0;        V1X = 0;
        V1Y = 0;        XC = 0;         YC = 0;

        for I = 1:NPE
           K  = (I-1)*NDF + 1;
           XC = XC + SFL(I)*ELXY(I,1);
           YC = YC + SFL(I)*ELXY(I,2);
           U1X = U1X + ELS(K)*GDSFL(1,I);
           U1Y = U1Y + ELS(K)*GDSFL(2,I);
           V1X = V1X + ELS(K+1)*GDSFL(1,I);
           V1Y = V1Y + ELS(K+1)*GDSFL(2,I);
        end
        
        XYC_ARRAY(COUNTER,1) = XC;
        XYC_ARRAY(COUNTER,2) = YC;

        %% Euler strain and Cauchy Stress tensor
        STRAIN = [U1X - 0.5*(U1X^2 + V1X^2) ; V1Y - 0.5*(U1Y^2 + V1Y^2) ; U1Y + V1X - (U1X*U1Y + V1X*V1Y)];
        CAUCHY_STRESS(COUNTER,:) = ([C(1,1) C(1,2) 0; C(1,2) C(2,2) 0; 0 0 C(6,6)]*STRAIN)';
        EULER_STRAIN(COUNTER,:) = STRAIN';
        
        ELXYC = zeros(size(ELXY));

        for I = 1:NPE
            L = (I-1)*NDF + 1;
            ELXYC(I,1) = ELXY(I,1) - ELS(L);
            ELXYC(I,2) = ELXY(I,2) - ELS(L+1);
        end

        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXYC);

        U1X = 0;        U1Y = 0;        V1X = 0;
        V1Y = 0;        X0 = 0;         Y0 = 0; 

        for I = 1:NPE
            L = (I-1)*NDF + 1;
            X0 = X0 + SFL(I)*ELXYC(I,1);
            Y0 = Y0 + SFL(I)*ELXYC(I,2);
            U1X = U1X + ELS(L)*GDSFL(1,I);
            U1Y = U1Y + ELS(L)*GDSFL(2,I);
            V1X = V1X + ELS(L+1)*GDSFL(1,I);
            V1Y = V1Y + ELS(L+1)*GDSFL(2,I);
        end
        
        XY0_ARRAY(COUNTER,1) = X0;
        XY0_ARRAY(COUNTER,2) = Y0;

        %% Green-Lagrange Strain and 2nd Piola Kirchhoff Stress tensor
        STRAIN = [U1X + 0.5*(U1X^2 + V1X^2) ; V1Y + 0.5*(U1Y^2 + V1Y^2) ; U1Y + V1X + (U1X*U1Y + V1X*V1Y)];
        PIOLA_STRESS(COUNTER,:) = ([C(1,1) C(1,2) 0; C(1,2) C(2,2) 0; 0 0 C(6,6)]*STRAIN)';
        GREEN_STRAIN(COUNTER,:) = STRAIN';

        COUNTER = COUNTER+1;
    end
end

%SS = [XY0_ARRAY PIOLA_STRESS GREEN_STRAIN zeros(size(XY0_ARRAY)) XYC_ARRAY CAUCHY_STRESS EULER_STRAIN];
SS = [CAUCHY_STRESS PIOLA_STRESS]*10^-5;

% SS.XY0 = XY0_ARRAY;
% SS.XYC = XYC_ARRAY;
% SS.PIOLA = PIOLA_STRESS;
% SS.GREEN = GREEN_STRAIN;
% SS.EULER = EULER_STRAIN;
% SS.CAUCHY = CAUCHY_STRESS;



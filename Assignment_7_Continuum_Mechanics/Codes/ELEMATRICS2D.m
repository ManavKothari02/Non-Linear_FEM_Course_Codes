function [ELK,ELF] = ELEMATRICS2D(NDF,NPE,ELXY,ELS,NGPF,C,thick,F,LFORM)


NN = NDF*NPE;
FX = 0;     FY = 0;
% Initialize ELK & ELF
ELF = zeros(NN,1);
ELK = zeros(NN);


[GAUSPT,GAUSWT] = Gaus_int(NGPF);

for NI = 1:NGPF
    for NJ = 1:NGPF
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ)*thick;

        X = 0;      Y = 0;      U = 0;
        V = 0;      DUX = 0;    DUY = 0;
        DVX = 0;    DVY = 0;
        
        for I = 1:NPE
            L = (I-1)*NDF + 1;
            X = X + SFL(I)*ELS(L);
            U = U + ELS(L)*SFL(I);
            DUX = DUX + GDSFL(1,I)*ELS(L);
            DUY = DUY + GDSFL(2,I)*ELS(L);
            DVX = DVX + GDSFL(1,I)*ELS(L+1);
            DVY = DVY + GDSFL(2,I)*ELS(L+1);
        end

        DUX2 = DUX^2;     DUY2 = DUY^2;
        DVX2 = DVX^2;     DVY2 = DVY^2;

        if LFORM == 1
            DUXP1 = 1 + DUX;
            DUXP2 = DUXP1^2;
            DVYP1 = 1 + DVY;
            DVYP2 = DVYP1^2;

            % Defining Green Strain and Second Piola Kirchhoff Stress
            STRAIN = zeros(3,1);

            STRAIN(1) = DUX + 0.5*(DUX2 + DVX2);
            STRAIN(2) = DVY + 0.5*(DUY2 + DVY2);
            STRAIN(3) = DUY + DVX + DUX*DUY + DVX*DVY;
            
            S11 = C(1,1)*STRAIN(1) + C(1,2)*STRAIN(2);
            S22 = C(1,2)*STRAIN(1) + C(2,2)*STRAIN(2);
            S12 = C(6,6)*STRAIN(3);

        else
            % Define Euler-Almansi Strain and Cauchy Stress Tensor
            STRAIN = zeros(3,1);

            STRAIN(1) = DUX - 0.5*(DUX2 + DVX2);
            STRAIN(2) = DVY - 0.5*(DUY2 + DVY2);
            STRAIN(3) = DUY + DVX - DUX*DUY - DVX*DVY;

            S11 = C(1,1)*STRAIN(1) + C(1,2)*STRAIN(2);
            S22 = C(1,2)*STRAIN(1) + C(2,2)*STRAIN(2);
            S12 = C(6,6)*STRAIN(3);

        end

        II = 1;
        
        for I = 1:NPE
            if LFORM == 1 % Total Lagrange
                ELF(II) = ELF(II) + (FX*F(1)*SFL(I) - DUXP1*GDSFL(1,I)*S11 - DUY*GDSFL(2,I)*S22 - (DUXP1*GDSFL(2,I) + DUY*GDSFL(1,I))*S12)*CONST;
                ELF(II+1) = ELF(II+1) + (FY*F(2)*SFL(I) - DVX*GDSFL(1,I)*S11 - DVYP1*GDSFL(2,I)*S22 - (DVYP1*GDSFL(1,I) + DVX*GDSFL(2,I))*S12)*CONST;
            else % Updated Lagrange
                ELF(II) = ELF(II) + (FX*F(1)*SFL(I) - GDSFL(1,I)*S11 - GDSFL(2,I)*S12)*CONST;
                ELF(II+1) = ELF(II+1) + (FY*F(2)*SFL(I) - GDSFL(1,I)*S12 - GDSFL(2,I)*S22)*CONST;
            end

            JJ = 1;

            for J = 1:NPE
                SIG = S11*GDSFL(1,I)*GDSFL(1,J) + S22*GDSFL(2,I)*GDSFL(2,J) + S12*(GDSFL(1,I)*GDSFL(2,J) + GDSFL(1,J)*GDSFL(2,I));

                if LFORM == 1 % Total Lagrange
                else % Updated Lagrange
                    ELK(II,JJ) = ELK(II,JJ) + (C(1,1)*GDSFL(1,I)*GDSFL(1,J) + C(6,6)*GDSFL(2,I)*GDSFL(2,J) + SIG)*CONST;
                    ELK(II+1,JJ+1) = ELK(II+1,JJ+1) + (C(6,6)*GDSFL(1,I)*GDSFL(1,J) + C(2,2)*GDSFL(2,I)*GDSFL(2,J) + SIG)*CONST;
                    ELK(II,JJ+1) = ELK(II,JJ+1) + (C(1,2)*GDSFL(1,I)*GDSFL(2,J) + C(6,6)*GDSFL(2,I)*GDSFL(1,J))*CONST;
                    ELK(II+1,JJ) = ELK(II+1,JJ) + (C(1,2)*GDSFL(2,I)*GDSFL(1,J) + C(6,6)*GDSFL(1,I)*GDSFL(2,J))*CONST;
                end
                JJ = NDF*J + 1;
            end
            II = NDF*I + 1;
        end
    end
end


















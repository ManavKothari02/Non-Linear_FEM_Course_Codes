function [ELK,ELF] = FLUIDMATRICS(NDF,NPE,NONLIN,ELXY,ELU,RHOAMU,NGPF,GAMA2)

RHO = RHOAMU(1);        MU = RHOAMU(2);

NN = NDF*NPE;
NGPR = NGPF-1;

% Initialize ELF, ELF and Tangent matrices

ELF = zeros(NN,1);
ELK = zeros(NN);
if NONLIN > 1
    TANG = zeros(NN);
end

[GAUSPT,GAUSWT] = Gaus_int(NGPF);

%% Full integration loop
for NI = 1:NGPF
    for NJ = 1:NGPF
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ);

        X = 0;      Y = 0;
        
        for I = 1:NPE
            X = X + ELXY(I,1)*SFL(I);
            Y = Y + ELXY(I,2)*SFL(I);
        end

        if NONLIN > 0 
            U = 0;      V = 0;
            DUX = 0;    DUY = 0;
            DVX = 0;    DVY = 0;

            for I = 1:NPE
                L = (I-1)*NDF + 1;
                U = ELU(L)*SFL(I) + U;
                V = ELU(L+1)*SFL(I) + V;
                DUX = ELU(L)*GDSFL(1,I) + DUX;
                DUY = ELU(L)*GDSFL(2,I) + DUY;
                DVX = ELU(L+1)*GDSFL(1,I) + DVX;
                DVY = ELU(L+1)*GDSFL(2,I) + DVY;
            end
        end

        II = 1;
        for I = 1:NPE
            JJ = 1;
            for J = 1:NPE
                S11 = GDSFL(1,I)*GDSFL(1,J)*CONST;
                S12 = GDSFL(1,I)*GDSFL(2,J)*CONST;
                S21 = GDSFL(2,I)*GDSFL(1,J)*CONST;
                S22 = GDSFL(2,I)*GDSFL(2,J)*CONST;

                ELK(II,JJ) = ELK(II,JJ) + MU*(2*S11 + S22);
                ELK(II+1,JJ) = ELK(II+1,JJ) + MU*S12;
                ELK(II,JJ+1) = ELK(II,JJ+1) + MU*S21;
                ELK(II+1,JJ+1) = ELK(II+1,JJ+1) + MU*(S11 + 2*S22);

                if NONLIN > 0
                    CNST = SFL(I)*(U*GDSFL(1,J) + V*GDSFL(2,J))*CONST;

                    ELK(II,JJ) = ELK(II,JJ) + RHO*CNST;
                    ELK(II+1,JJ+1) = ELK(II+1,JJ+1) + RHO*CNST;

                    if NONLIN > 1
                        TANG(II,JJ) = RHO*DUX*S00;
                        TANG(II,JJ+1) = RHO*DUY*S00;
                        TANG(II+1,JJ) = RHO*DVX*S00;
                        TANG(II+1,JJ+1) = RHO*DVY*S00;
                    end
                end
                JJ = NDF*J + 1;
            end
            II = NDF*I + 1;
        end
    end
end % End of full integration loop

%% Reduced integration for Penalty Terms
[GAUSPT,GAUSWT] = Gaus_int(NGPR);

for NI = 1:NGPR
    for NJ = 1:NGPR
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ);

        II = 1;
        for I = 1:NPE
            JJ = 1;
            for J = 1:NPE
                S11 = GDSFL(1,I)*GDSFL(1,J)*CONST;
                S12 = GDSFL(1,I)*GDSFL(2,J)*CONST;
                S21 = GDSFL(2,I)*GDSFL(1,J)*CONST;
                S22 = GDSFL(2,I)*GDSFL(2,J)*CONST;

                ELK(II,JJ) = ELK(II,JJ) + GAMA2*S11;
                ELK(II+1,JJ) = ELK(II+1,JJ) + GAMA2*S21;
                ELK(II,JJ+1) = ELK(II,JJ+1) + GAMA2*S12;
                ELK(II+1,JJ+1) = ELK(II+1,JJ+1) + GAMA2*S22;

                JJ = NDF*J + 1;
            end
            II = I*NDF + 1;
        end
    end
end % End of Reduced Integration Loop

%% Calculation of Tangent and Residual Matrices

if NONLIN > 1
    for I = 1:NN
        for J = 1:NN
            ELF(I) = ELF(I) - ELK(I,J)*ELU(J);

            ELK(I,J) = ELK(I,J) + TANG(I,J);
        end
    end
end
end












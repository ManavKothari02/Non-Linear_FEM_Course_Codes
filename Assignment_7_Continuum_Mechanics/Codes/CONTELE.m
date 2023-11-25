function [ELK,ELF] = CONTELE(NDF,NPE,ELXY,ELS,NGPF,C,thick,F)

NN = NDF*NPE;

% Initialize ELK & ELF

K111 = zeros(NPE);      K121 = zeros(NPE);
K112 = zeros(NPE);      K211 = zeros(NPE);
K221 = zeros(NPE);      K222 = zeros(NPE);

FORCE211 = zeros(NPE,1);    FORCE212 = FORCE211;
FORCE111 = FORCE212;        FORCE112 = FORCE111;

[GAUSPT,GAUSWT] = Gaus_int(NGPF);


%% Full integration loop

for NI = 1:NGPF
    for NJ = 1:NGPF
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ);

        X = 0;      Y = 0;      U = 0;
        V = 0;      DUX = 0;    DUY = 0;
        DVX = 0;    DVY = 0;

        for I = 1:NPE
            L = (I-1)*NDF + 1;
            X = X + ELXY(I,1)*SFL(I);
            Y = Y + ELXY(I,2)*SFL(I);
            U = U + ELS(L)*SFL(I);
            V = V + ELS(L+1)*SFL(I);
            DUX = DUX + ELS(L)*GDSFL(1,I);
            DUY = DUY + ELS(L)*GDSFL(2,I);
            DVX = DVX + ELS(L+1)*GDSFL(1,I);
            DVY = DVY + ELS(L+1)*GDSFL(2,I);
        end
        
        % Defining Green Strain & 2nd Piola - Kirchhoff stress components
        STRAIN = [DUX - 0.5*(DUX^2) - 0.5*(DVX^2) ; DVY - 0.5*(DUY^2) - 0.5*(DVY^2) ; DUY + DVX - (DUX*DUY) - (DVX*DVY)]; 
        SIGMA = [C(1,1) C(1,2) 0; C(1,2) C(2,2) 0; 0 0 C(6,6)]*STRAIN;
        
        for I = 1:NPE
            % Calculation of force vector
            FORCE211(I) = FORCE211(I) + thick*F(1)*SFL(I)*CONST;
            FORCE212(I) = FORCE212(I) + thick*F(2)*SFL(I)*CONST;
            FORCE111(I) = FORCE111(I) + thick*(SIGMA(1)*GDSFL(1,I) + SIGMA(3)*GDSFL(2,I))*CONST;
            FORCE112(I) = FORCE112(I) + thick*(SIGMA(3)*GDSFL(1,I) + SIGMA(2)*GDSFL(2,I))*CONST;

            for J = 1:NPE
                S11 = GDSFL(1,I)*GDSFL(1,J)*CONST;
                S12 = GDSFL(1,I)*GDSFL(2,J)*CONST;
                S21 = GDSFL(2,I)*GDSFL(1,J)*CONST;
                S22 = GDSFL(2,I)*GDSFL(2,J)*CONST;

                K111(I,J) = K111(I,J) + thick*(C(1,1)*S11 + C(6,6)*S22);
                K121(I,J) = K121(I,J) + thick*(C(1,2)*S12 + C(6,6)*S21);
                K211(I,J) = K211(I,J) + thick*(C(1,2)*S21 + C(6,6)*S12);
                K221(I,J) = K221(I,J) + thick*(C(6,6)*S11 + C(2,2)*S22);
                K112(I,J) = K112(I,J) + thick*(SIGMA(1)*S11 + SIGMA(3)*(S21+S12) + SIGMA(2)*S22);
                K222(I,J) = K112(I,J);
                
            end
        end
    end
end % End of Full-integration loops

ELK = [K111+K112 K121 ; K211 K221+K222];
ELF = [FORCE211 - FORCE111 ; FORCE212 - FORCE112];

%% Calculation of Tangent and Residual Matrices
for I = 1:NN
    for J = 1:NN
        ELF(I) = ELF(I) - ELK(I,J)*ELS(J);
        %ELK(I,J) = ELK(I,J) + TANG(I,J);
    end
end
end
















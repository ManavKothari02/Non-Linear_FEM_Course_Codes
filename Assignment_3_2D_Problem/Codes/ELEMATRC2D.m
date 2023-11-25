function [ELK,ELF] = ELEMATRC2D(NGPF,NPE,NONLIN,ELXY,ELU,PDECOEFF)

[GAUSPT,GAUSWT] = Gaus_int(NGPF);

A11 = PDECOEFF.A11;
A22 = PDECOEFF.A22;
A00 = PDECOEFF.A00;
F = PDECOEFF.F;

A10 = A11(1);   A1X = A11(2);   A1Y = A11(3);
A1U = A11(4);   A1UX = A11(5);  A1UY = A11(6);

A20 = A22(1);   A2X = A22(2);   A2Y = A22(3);
A2U = A22(4);   A2UX = A22(5);  A2UY = A22(6);

A0X = A00(2);   A0Y = A00(3);   A00 = A00(1);

F0 = F(1);      FX = F(2);      FY = F(3);

% Initialize Matrices

ELK = zeros(NPE);
ELF = zeros(NPE,1);

if NONLIN > 1 % ie, Newton's iteration
    TANG = zeros(NPE);
end

%% Loop for numerical integration
for NI = 1:NGPF
    for NJ = 1:NGPF

        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ);

        X = 0; Y = 0; U = 0; DUX = 0; DUY = 0;
        
        for I = 1:NPE
            if NONLIN > 0
                U = U + ELU(I)*SFL(I);
                DUX = DUX + ELU(I)*GDSFL(1,I);
                DUY = DUY + ELU(I)*GDSFL(2,I);
            end
            X = X + ELXY(I,1)*SFL(I);
            Y = Y + ELXY(I,2)*SFL(I);
        end

        % Define coefficients of the Differential Equation
        
        FXY = F0 + FX*X + FY*Y;
        A11 = A10 + A1X*X + A1Y*Y;
        A22 = A20 + A2X*X + A2Y*Y;

        if NONLIN > 0 
            AXX = A11 + A1U*U + A1UX*DUX + A1UY*DUY;
            AYY = A22 + A2U*U + A2UX*DUX + A2UY*DUY;
        end

        % Defining ELF and ELF
        for I = 1:NPE
            ELF(I) = ELF(I) + SFL(I)*CONST*FXY;
            for J = 1:NPE
                S00 = SFL(I)*SFL(J)*CONST;
                S11 = GDSFL(1,I)*GDSFL(1,J)*CONST;
                S22 = GDSFL(2,I)*GDSFL(2,J)*CONST;
                ELK(I,J) = ELK(I,J) + A00*S00 + AXX*S11 + AYY*S22;
            end
        end
        if NONLIN > 1
            S10 = GDSFL(1,I)*SFL(J)*CONST;
            S20 = GDSFL(2,I)*SFL(J)*CONST;
            S12 = GDSFL(1,I)*GDSFL(2,J)*CONST;
            S21 = GDSFL(2,I)*GDSFL(1,J)*CONST;
            
            TANG(I,J) = TANG(I,J) + DUX*(A1U*S10+A1UX*S11+A1UY*S12) + DUY*(A2U*S20+A2UX*S21+A2UY*S22);
        end
    end
end

% Compute Residual Vector and tangent matrix

if NONLIN > 1
    for I = 1:NPE
        for J = 1:NPE
            ELF(I) = ELF(I) - ELK(I,J)*ELU(J);

            ELK(I,J) = ELK(I,J) + TANG(I,J);
        end
    end
end

end
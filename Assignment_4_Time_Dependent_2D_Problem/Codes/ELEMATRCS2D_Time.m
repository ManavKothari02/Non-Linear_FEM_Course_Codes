function [ELK,ELM,ELF] = ELEMATRCS2D_Time(NGPF,NPE,NDF,NONLIN,ELXY,ELU,ELU0,ELU1,ELU2,PDECOEFF,ITEM,alfa,GAMA,DT,A_Coeff)
DT2 = DT^2;
NN = NPE*NDF;

A1 = A_Coeff(1);    A2 = A_Coeff(2);    A3 = A_Coeff(3);    A4 = A_Coeff(4);
A5 = A_Coeff(5);    A6 = A_Coeff(6);    A7 = A_Coeff(7);    A8 = A_Coeff(8);

A11 = PDECOEFF.A11;     A22 = PDECOEFF.A22;
C = PDECOEFF.C;         A00 = PDECOEFF.A00;
f = PDECOEFF.F;

A10 = A11(1);   A1x = A11(2);   A1y = A11(3);
A1u = A11(4);   A1ux = A11(5);  A1uy = A11(6);

A20 = A22(1);   A2x = A22(2);   A2y = A22(3);
A2u = A22(4);   A2ux = A22(5);  A2uy = A22(6);

A0x = A00(2);   A0y = A00(3);   A0u = A00(4);
A0ux = A00(5);  A0uy = A00(6);  A00 = A00(1);

C0 = C(1);      Cx = C(2);      Cy = C(3);

F0 = f(1);      Fx = f(2);      Fy = f(3);

% Initialize ELK, ELF and ELM
ELK = zeros(NPE);
ELM = zeros(NPE);
ELF = zeros(NPE,1);
ELK0 = zeros(NPE);

if NONLIN > 1 
    TANG = zeros(NPE);
end

%% Loops for numerical integration

[GAUSPT,GAUSWT] = Gaus_int(NGPF);

for NI = 1:NGPF
    for NJ = 1:NGPF

        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ);

        x = 0;  y = 0;  u = 0;  dux = 0;    duy = 0;

        for I = 1:NPE
            if NONLIN > 0
                u = u + ELU(I)*SFL(I);
                dux = dux + ELU(I)*GDSFL(1,I);
                duy = duy + ELU(I)*GDSFL(2,I);
            end
            x = x + ELXY(I,1)*SFL(I);
            y = y + ELXY(I,2)*SFL(I);
        end

        Fxy = F0 + Fx*x + Fy*y;
        A11 = A10 + A1x*x + A1y*y;
        A22 = A20 + A2x*x + A2y*y;
        A00 = A00 + A0x*x + A0y*y;

        if NONLIN > 0
            Axx = A11 + A1u*u + A1ux*dux + A1uy*duy;
            Ayy = A22 + A2u*u + A2ux*dux + A2uy*duy;
            A00 = A00 + A0u*u + A0ux*dux + A0uy*duy;
        end
        
        if ITEM > 0
            Cxy = C0 + Cx*x + Cy*y;
        end
        
        if ITEM > 0
            UP = 0;     UPx = 0;    UPy = 0;
            for I = 1:NPE
                UP = UP + ELU0(I)*SFL(I);
                UPx = UPx + ELU0(I)*GDSFL(1,I);
                UPy = UPy + ELU0(I)*GDSFL(2,I);
            end
            APxx = A11 + A1u*UP + A1ux*UPx + A1y*UPy;
            APyy = A22 + A2u*UP + A2ux*UPx + A2uy*UPy;
        end

        % Defining ELF and ELK
        for I = 1:NPE
            ELF(I) = ELF(I) + SFL(I)*CONST*Fxy;
            for J = 1:NPE
                S00 = SFL(I)*SFL(J)*CONST;
                S11 = GDSFL(1,I)*GDSFL(1,J)*CONST;
                S22 = GDSFL(2,I)*GDSFL(2,J)*CONST;

                ELK(I,J) = ELK(I,J) + Axx*S11 + Ayy*S22 + A00*S00;
                
                if NONLIN > 1
                    S10 = GDSFL(1,I)*SFL(J)*CONST;
                    S20 = GDSFL(2,I)*SFL(J)*CONST;
                    S12 = GDSFL(1,I)*GDSFL(2,I)*CONST;
                    S21 = GDSFL(2,I)*GDSFL(1,I)*CONST;
                    TANG(I,J) = TANG(I,J) + dux*(A1u*S10 + A1ux*S11 + A1uy*S12) + duy*(A2u*S20 + A2ux*S21 + A2uy*S22);
                end

                if (ITEM > 0)
                    ELM(I,J) = ELM(I,J) + Cxy*S00;
                    if NONLIN > 0
                        ELK0(I,J) = ELK0(I,J) + APxx*S11 + APyy*S22 + A00*S00;
                    end
                end
            end
        end
    end
end

if ITEM == 1 % Parabolic Equation
    for I = 1:NN
        SUM = 0;
        for J = 1:NN
            if NONLIN > 0
                SUM = SUM + (ELM(I,J) - A2*ELK0(I,J))*ELU0(J);
            else
                SUM = SUM + (ELM(I,J)- A2*ELK(I,J)) * ELU0(J);
            end
            ELK(I,J) = ELM(I,J) + A1*ELK(I,J);
        end
            ELF(I) = (A1+A2)*ELF(I) + SUM;
    end
elseif ITEM > 1 % Hyperbolic Equation
    for I = 1:NN
        SUM = 0;
        for J = 1:NN
            SUM = SUM + ELM(I,J)*(A3*ELU0(J) + A4*ELU1(J) + A5*ELU2(J));
            ELK(I,J) = ELK(I,J) + A3*ELM(I,J);
        end
        ELF(I) = ELF(I) + SUM;
    end
end

if NONLIN > 1
    ELF = ELF - ELK*ELU;
    ELK = ELK + TANG;
end


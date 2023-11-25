function [ELK,ELF] = FSDTPLTS(NGPF,NGPR,NGPS,NPE,NDF,NONLIN,ELXY,ELU,q0,E,v,G,Ks,h)

NN = NPE*NDF; % Number of DOF per node

qx = 0;         qy = 0;

E1 = E(1);      E2 = E(2);
v12 = v(1);     v21 = v(2);
G12 = G(1);     G13 = G(2);     G23 = G(3);

A11 = E1*h/(1-v12*v21);     A12 = v21*A11;
A55 = Ks*G13*h;             A66 = G12*h;
A44 = Ks*G23*h;             A22 = A11*(E2/E1);
D11 = A11*(h^2/12);         D12 = D11*v21;
D22 = D11*(E2/E1);          D66 = G12*(h^3/12);


ELK = zeros(NN);
ELF = zeros(NN,1);

if NONLIN > 1
    TANG = zeros(NN);
end

[GAUSPT,GAUSWT] = Gaus_int(NGPF);

%% Loop for numerical integration

for NI = 1:NGPF
    for NJ = 1:NGPF
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ);

        X = 0; Y = 0;

        for I = 1:NPE
            X = X + ELXY(I,1)*SFL(I);
            Y = Y + ELXY(I,2)*SFL(I);
        end

        QXY = q0 + qx*X + qy*Y;

        II = 1;
        for I = 1:NPE
            II1 = II+1;
            II2 = II+2;
            II3 = II+3;
            II4 = II+4;

            ELF(II2) = ELF(II2) + CONST*SFL(I)*QXY;
                
            JJ = 1;
            for J = 1:NPE
                JJ1 = JJ+1;
                JJ2 = JJ+2;
                JJ3 = JJ+3;
                JJ4 = JJ+4;

                SXX = GDSFL(1,I)*GDSFL(1,J);
                SXY = GDSFL(1,I)*GDSFL(2,J);
                SYX = GDSFL(2,I)*GDSFL(1,J);
                SYY = GDSFL(2,I)*GDSFL(2,J);

                % Stiffness coeff for inplane DOF
                ELK(II,JJ) = ELK(II,JJ) + CONST*(A11*SXX + A66*SYY);
                ELK(II,JJ1) = ELK(II,JJ1) + CONST*(A12*SXY + A66*SYX);
                ELK(II1,JJ) = ELK(II1,JJ) + CONST*(A12*SYX + A66*SXY);
                ELK(II1,JJ1) = ELK(II1,JJ1) + CONST*(A66*SXX + A22*SYY);

                % Bending Stiffness
                ELK(II3,JJ3) = ELK(II3,JJ3) + CONST*(D11*SXX + D66*SYY);
                ELK(II3,JJ4) = ELK(II3,JJ4) + CONST*(D12*SXY + D66*SYX);
                ELK(II4,JJ3) = ELK(II4,JJ3) + CONST*(D12*SYX + D66*SXY);
                ELK(II4,JJ4) = ELK(II4,JJ4) + CONST*(D66*SXX + D22*SYY);

                JJ = NDF*J + 1;
            end
            II = NDF*I + 1;
        end
    end
end

%% Reduced Integration to evaluate transverse shear stiffness

[GAUSPT,GAUSWT] = Gaus_int(NGPS);

for NI = 1:NGPS
    for NJ = 1:NGPS
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        CONST = JAC*GAUSWT(NI)*GAUSWT(NJ);
        II = 1;
        for I = 1:NPE
            II1 = II+1;
            II2 = II+2;
            II3 = II+3;
            II4 = II+4;

            JJ = 1;
            for J = 1:NPE
                JJ1 = JJ+1;
                JJ2 = JJ+2;
                JJ3 = JJ+3;
                JJ4 = JJ+4;

                S00 = SFL(I)*SFL(J);
                SX0 = GDSFL(1,I)*SFL(J);
                S0X = SFL(I)*GDSFL(1,J);
                SY0 = GDSFL(2,I)*SFL(J);
                S0Y = SFL(I)*GDSFL(2,J);
                SXX = GDSFL(1,I)*GDSFL(1,J);
                SXY = GDSFL(1,I)*GDSFL(2,J);
                SYX = GDSFL(2,I)*GDSFL(1,J);
                SYY = GDSFL(2,I)*GDSFL(2,J);

                ELK(II2,JJ2) = ELK(II2,JJ2) + CONST*(A55*SXX + A44*SYY);
                ELK(II2,JJ3) = ELK(II2,JJ3) + CONST*A55*SX0;
                ELK(II3,JJ2) = ELK(II3,JJ2) + CONST*A55*S0X;
                ELK(II2,JJ4) = ELK(II2,JJ4) + CONST*A44*SY0;
                ELK(II4,JJ2) = ELK(II4,JJ2) + CONST*A44*S0Y;
                ELK(II3,JJ3) = ELK(II3,JJ3) + CONST*A55*S00;
                ELK(II4,JJ4) = ELK(II4,JJ4) + CONST*A44*S00;

                JJ = NDF*J + 1;
            end
            II = NDF*I + 1;
        end
    end
end

%% Reduced integration for nonlinear stiffness

[GAUSPT,GAUSWT] = Gaus_int(NGPR);
if NONLIN > 0
    for NI = 1:NGPR
        for NJ = 1:NGPR
            [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
            CONST = JAC*GAUSWT(NI)*GAUSWT(NJ); 

            DUX = 0;        DUY = 0;
            DVX = 0;        DVY = 0;
            DWX = 0;        DWY = 0;
            DPXX = 0;       DPXY = 0;
            DPYX = 0;       DPYY = 0;
            
            for I = 1:NPE
                DUX = DUX + GDSFL(1,I)*ELU(I*NDF-4);
                DUY = DUY + GDSFL(2,I)*ELU(I*NDF-4);
                DVX = DVX + GDSFL(1,I)*ELU(I*NDF-3);
                DVY = DVY + GDSFL(2,I)*ELU(I*NDF-3);
                DWX = DWX + GDSFL(1,I)*ELU(I*NDF-2);
                DWY = DWY + GDSFL(2,I)*ELU(I*NDF-2);
                DPXX = DPXX + GDSFL(1,I)*ELU(I*NDF-1);
                DPXY = DPXY + GDSFL(2,I)*ELU(I*NDF-1);
                DPYX = DPYX + GDSFL(1,I)*ELU(I*NDF);
                DPYY = DPYY + GDSFL(2,I)*ELU(I*NDF);
            end

            DWX2 = DWX^2;
            DWY2 = DWY^2;
            DWXY = DWX*DWY;
            EXX_0 = DUX + 0.5*DWX2;
            EYY_0 = DVY + 0.5*DWY2;
            GXY_0 = DUY + DVX + DWXY;
            NXX = A11*EXX_0 + A12*EYY_0;
            NYY = A12*EXX_0 + A22*EYY_0;
            NXY = A66*GXY_0;
            II = 1;

            for I = 1:NPE
                II1 = II+1;
                II2 = II+2;
                II3 = II+3;
                II4 = II+4;
                JJ = 1;

                for J = 1:NPE
                    JJ1 = JJ+1;
                    JJ2 = JJ+2;
                    JJ3 = JJ+3;
                    JJ4 = JJ+4;

                    SXX = GDSFL(1,I)*GDSFL(1,J);
                    SXY = GDSFL(1,I)*GDSFL(2,J);
                    SYX = GDSFL(2,I)*GDSFL(1,J);
                    SYY = GDSFL(2,I)*GDSFL(2,J);

                    ELK(II,JJ2)=ELK(II,JJ2) + 0.5*CONST*(DWX*A11*SXX + DWY*A12*SXY + A66*(DWX*SYY + DWY*SYX));
                    ELK(II1,JJ2)=ELK(II1,JJ2) + 0.5*CONST*(DWX*A12*SYX + DWY*A22*SYY + A66*(DWX*SXY + DWY*SXX));
                    ELK(II2,JJ)=ELK(II2,JJ) + CONST*(DWX*A11*SXX + DWY*A12*SYX + A66*(DWX*SYY + DWY*SXY));
                    ELK(II2,JJ1)=ELK(II2,JJ1) + CONST*(DWX*A12*SXY + DWY*A22*SYY + A66*(DWX*SYX + DWY*SXX));
                    ELK(II2,JJ2)=ELK(II2,JJ2) + 0.5*CONST*(A11*SXX*DWX2 + A66*SXX*DWY2 + A66*SYY*DWX2 + A22*DWY2*SYY + (A12 + A66)*DWXY*(SXY + SYX));
                    
                    if NONLIN > 1
                        TANG(II,JJ2)=TANG(II,JJ2) + 0.5*CONST*(A11*DWX*SXX + A12*DWY*SXY + A66*(DWX*SYY + DWY*SYX));
                        %TANG(II2,JJ)=TANG(II2,JJ) + ELK(II2,JJ);
                        TANG(II1,JJ2)=TANG(II1,JJ2) + 0.5*CONST*(A12*DWX*SYX + A22*DWY*SYY + A66*(DWY*SXX + DWX*SXY));
                        %TANG(II2,JJ1)=TANG(II2,JJ1) + ELK(II2,JJ1);
                        TANG(II2,JJ2)=TANG(II2,JJ2) + CONST*(NXX*SXX + NYY*SYY + NXY*SXY + NXY*SYX + 0.5*((A12+A66)*DWXY*(SXY+SYX) ...
                                        +SXX*(A11*DWX2+A66*DWY2) + SYY*(A66*DWX2 + A22*DWY2)));
                    end
                    JJ = J*NDF + 1;
                end
                II = I*NDF+1;
            end
        end
    end
end

%% Compute Residual Vector and Tangent Matrix

if NONLIN > 1
    for I = 1:NN
        for J = 1:NN
            ELF(I) = ELF(I) - ELK(I,J)*ELU(J);

            ELK(I,J) = ELK(I,J) + TANG(I,J);
        end
    end
end
end



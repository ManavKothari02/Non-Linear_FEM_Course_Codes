function [ELK,ELF] = ELMATRCS1D(model, IEL, NONLIN, ELX, ELU, NGP, LGP, E, A, I, G, Ks, f, F, Q)

[GAUSPT,GAUSWT] = Gaus_int(NGP);

NPE = length(ELX);

ELF1 = zeros(NPE,1);
ELF2 = zeros(NPE,1);
ELF3 = zeros(NPE,1);

K11 = zeros(NPE);
K12 = K11;      K13 = K11;
K21 = K11;      K22 = K11;
K23 = K11;      K31 = K11;
K32 = K11;      K33 = K11;

% Initializing Tangent Matrix
TAN12 = zeros(NPE);
TAN13 = TAN12;      TAN22 = TAN12;
TAN23 = TAN12;      TAN32 = TAN12;
TAN33 = TAN12;

he = ELX(end) - ELX(1);
NDF = 3;

for NI = 1:NGP
    [SFL,GDSFL,SFH,GDSFH,GDDSFH,GJ] = interpolation_function(ELX,IEL,model,NPE,GAUSPT(NI));
    x = ELX(1) + 0.5*(1+GAUSPT(NI))*he;
    Const = GJ*GAUSWT(NI);

    % Defining Axx, Dxx, Sxz, Fx, Qx
    Axx = E*A;
    Dxx = E*I;
    Sxz = G*A*Ks;

    Fx = F(1) + F(2)*x + F(3)*x^2;
    Qx = Q(1) + Q(2)*x + Q(3)*x^2;

    %% Euler Bernoulli Beam Theory
    if model == 2
        for i = 1:NPE
            count1 = 2*i - 1;
            ELF1(i) = ELF1(i) + f*Fx*SFL(i)*Const;
            ELF2(i) = ELF2(i) + f*Qx*SFH(count1)*Const;
            ELF3(i) = ELF3(i) + f*Qx*SFH(count1+1)*Const;

            for j = 1:NPE
                count2 = 2*j-1;
                S11 = GDSFL(i)*GDSFL(j)*Const;
                H22 = GDDSFH(count1)*GDDSFH(count2)*Const;
                H23 = GDDSFH(count1)*GDDSFH(count2+1)*Const;
                H32 = GDDSFH(count1+1)*GDDSFH(count2)*Const;
                H33 = GDDSFH(count1+1)*GDDSFH(count2+1)*Const;

                K11(i,j) = K11(i,j) + Axx*S11;
                K22(i,j) = K22(i,j) + Dxx*H22;
                K23(i,j) = K23(i,j) + Dxx*H23;
                K32(i,j) = K32(i,j) + Dxx*H32;
                K33(i,j) = K33(i,j) + Dxx*H33;
            end
        end
    %% Timoshenko Beam Theory
    elseif model == 3
        for i = 1:NPE
            ELF1(i) = ELF1(i) + f*Fx*SFL(i)*Const;
            ELF2(i) = ELF2(i) + f*Qx*SFL(i)*Const;
            for j = 1:NPE
                S11 = GDSFL(i)*GDSFL(j)*Const;
                K11(i,j) = K11(i,j) + Axx*S11;
                K33(i,j) = K33(i,j) + Dxx*S11;
            end
        end
    end
end

% Reduced integration
U = 0;
dU = 0;
dw = 0;

[GAUSPT,GAUSWT] = Gaus_int(LGP);

for k = 1:(LGP)
    [SFL,GDSFL,SFH,GDSFH,GDDSFH,GJ] = interpolation_function(ELX,IEL,model,NPE,GAUSPT(k));
    % x = ELX(1) + 0.5*(1+GAUSPT(k))*he;
    Const = GJ*GAUSWT(k);
    for i = 1:NPE
        L = i*NDF-1;
        K = 2*i-1;
        dU = dU + GDSFL(i)*ELU(L-1);
        if model == 2
            dw = dw + GDSFH(K)*ELU(L) + GDSFH(K+1)*ELU(L+1);
        elseif model == 3
            dw = dw + GDSFL(i)*ELU(L);
        end
    end

    AXW = Axx*dw;
    AXWH = 0.5*AXW;
    AN0 = AXWH*dw;
    AN1 = Axx*(dU + dw*dw);

    if model == 2 % Euler Bernoulli Non-Linear Terms
        for i = 1:NPE
            count1 = 2*i-1;
            for j = 1:NPE
                count2 = 2*j-1;
                K12(i,j) = K12(i,j) + AXWH*GDSFL(i)*GDSFH(count2)*Const;
                K13(i,j) = K13(i,j) + AXWH*GDSFL(i)*GDSFH(count2+1)*Const;
                K21(i,j) = K21(i,j) + AXW*GDSFL(j)*GDSFH(count1)*Const;
                K31(i,j) = K31(i,j) + AXW*GDSFL(j)*GDSFH(count1+1)*Const;
                K22(i,j) = K22(i,j) + AN0*GDSFH(count1)*GDSFH(count2)*Const;
                K23(i,j) = K23(i,j) + AN0*GDSFH(count1)*GDSFH(count2+1)*Const;
                K32(i,j) = K32(i,j) + AN0*GDSFH(count1+1)*GDSFH(count2)*Const;
                K33(i,j) = K33(i,j) + AN0*GDSFH(count1+1)*GDSFH(count2+1)*Const;

                if NONLIN == 2
                    TAN12(i,j) = TAN12(i,j) + AXWH*GDSFL(i)*GDSFH(count2)*Const; 
                    TAN13(i,j) = TAN13(i,j) + AXWH*GDSFL(i)*GDSFH(count2+1)*Const;
                    TAN22(i,j) = TAN22(i,j) + AN1*GDSFH(count1)*GDSFH(count2)*Const;
                    TAN23(i,j) = TAN23(i,j) + AN1*GDSFH(count1)*GDSFH(count2+1)*Const;
                    TAN32(i,j) = TAN32(i,j) + AN1*GDSFH(count1+1)*GDSFH(count2)*Const;
                    TAN33(i,j) = TAN33(i,j) + AN1*GDSFH(count1+1)*GDSFH(count2+1)*Const;
                end
            end
        end
    elseif model == 3 % Timoshenko Beam Theory
        for i = 1:NPE
            for j = 1:NPE
                S11 = GDSFL(i)*GDSFL(j)*Const;
                S01 = SFL(i)*GDSFL(j)*Const;
                S10 = GDSFL(i)*SFL(j)*Const;
                S00 = SFL(i)*SFL(j)*Const;

                K22(i,j) = K22(i,j) + Sxz*S11;
                K23(i,j) = K23(i,j) + Sxz*S10;
                K32(i,j) = K32(i,j) + Sxz*S01;
                K33(i,j) = K33(i,j) + Sxz*S00;

                if NONLIN > 0
                    K12(i,j) = K12(i,j) + 0.5*AXW*S11;
                    K21(i,j) = K21(i,j) + AXW*S11;
                    K22(i,j) = K22(i,j) + AN0*S11;
                    if NONLIN == 2
                        TAN12(i,j) = TAN12(i,j) + 0.5*AXW*S11;
                        TAN22(i,j) = TAN22(i,j) + AN1*S11;
                    end
                end
            end
        end
    end
end

%% Rearranging Element Coefficient
NET = NPE*NDF;

ELK = zeros(NET);
ELF = zeros(NET,1);
if model > 1
    Flag1 = 1;
    for i = 1:NPE
        ELF(Flag1) = ELF1(i);
        ELF(Flag1+1) = ELF2(i);
        ELF(Flag1+2) = ELF3(i);

        Flag2 = 1;
        for j = 1:NPE
            ELK(Flag1,Flag2) = K11(i,j);
            ELK(Flag1,Flag2+1) = K12(i,j);
            ELK(Flag1,Flag2+2) = K13(i,j);
            ELK(Flag1+1,Flag2) = K21(i,j);
            ELK(Flag1+2,Flag2) = K31(i,j);
            ELK(Flag1+1,Flag2+1) = K22(i,j);
            ELK(Flag1+1,Flag2+2) = K23(i,j);
            ELK(Flag1+2,Flag2+1) = K32(i,j);
            ELK(Flag1+2,Flag2+2) = K33(i,j);
            Flag2 = NDF*j+1;
        end
        Flag1 = NDF*i+1;
    end
end

%% Calculating Residual vector for Newton's Iteration

if NONLIN == 2
    ELF = ELF - ELK*ELU;

    II = 1;
    
    for i = 1:NPE
        JJ = 1;
        for j = 1:NPE
            ELK(II,JJ+1) = ELK(II,JJ+1) + TAN12(i,j);
            ELK(II+1,JJ+1) = ELK(II+1,JJ+1) + TAN22(i,j);
            if model>=2
                ELK(II,JJ+2) = ELK(II,JJ+2) + TAN13(i,j);
                ELK(II+1,JJ+2) = ELK(II+1,JJ+2) + TAN23(i,j);
                ELK(II+2,JJ+1) = ELK(II+2,JJ+1) + TAN32(i,j);
                ELK(II+2,JJ+2) = ELK(II+2,JJ+2) + TAN33(i,j);
                JJ = NDF*j+1;
            end
        end
        II = NDF*i+1;
    end
end
















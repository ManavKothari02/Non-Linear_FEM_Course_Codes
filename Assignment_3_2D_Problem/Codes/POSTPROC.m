function POSTPROC(ELXY,NDF,NPE,ELU,NGPF,NGPR,PDECOEFF,iter,convergence)
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

[GAUSPT,GAUSWT] = Gaus_int(NGPR);

for NI = 1:NGPR
    for NJ = 1:NGPR
        XI = GAUSPT(NI);
        ETA = GAUSPT(NJ);
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,XI,ETA,ELXY);
        XC = 0; YC = 0;

        for I = 1:NPE
            XC = XC + ELXY(I,1)*SFL(I);
            YC = YC + ELXY(I,2)*SFL(I);
        end

        DUX = 0; DUY = 0;

        for J = 1:NPE
            DUX = DUX + ELU(J)*GDSFL(1,J);
            DUY = DUY + ELU(J)*GDSFL(2,J);
        end

        FluxX = -(A10 + A1X*XC + A1Y*YC)*DUX;
        FluxY = -(A20 + A2X*XC + A2Y*YC)*DUY;
        NetFlux = sqrt(FluxX^2 + FluxY^2);

        if FluxX == 0
            if FluxY < 0
                ANGLE = -90;
            else
                ANGLE = 90;
            end
        else
            ANGLE = atan2(FluxX,FluxY)*180/pi;
        end
        if convergence == 1
            iter = iter-1
            XC
            YC
            FluxX
            FluxY
            NetFlux
            ANGLE
        end
        
        
    end
end




end

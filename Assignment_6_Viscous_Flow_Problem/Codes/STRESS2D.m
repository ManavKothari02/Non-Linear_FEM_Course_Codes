function [XMAT,YMAT,SXMAT,SYMAT,SXYMAT,PRSMAT] = STRESS2D(ELXY,NPE,ELU,NGPR,GAMA2,MU)

NDF = 2;

[GAUSPT,GAUSWT] = Gaus_int(NGPR);

for NI = 1:NGPR
    for NJ = 1:NGPR
        [SFL,GDSFL,JAC] = INTERPLN2D(NPE,GAUSPT(NI),GAUSPT(NJ),ELXY);
        X = 0;      Y = 0;      U = 0;      V = 0;
        UX = 0;     UY = 0;     VX = 0;     VY = 0;

        for I = 1:NPE
            X = X + SFL(I)*ELXY(I,1);
            Y = Y + SFL(I)*ELXY(I,2);

            L = (I-1)*NDF + 1;

            UX = UX + ELU(L)*GDSFL(1,I);
            UY = UY + ELU(L)*GDSFL(2,I);
            VX = VX + ELU(L+1)*GDSFL(1,I);
            VY = VY + ELU(L+1)*GDSFL(2,I);
        end

        DIV = UX + VY;
        PRS = -GAMA2*DIV;
        SX = 2*MU*UX - PRS;
        SY = 2*MU*VY - PRS;
        SXY = MU*(UY + VX);
        
        XMAT = [];
        YMAT = [];
        PRSMAT = [];
        SXMAT = [];
        SYMAT = [];
        SXYMAT = [];

        XMAT(length(XMAT)+1) = X;
        YMAT(length(YMAT)+1) = Y;
        PRSMAT(length(PRSMAT)+1) = PRS;
        SXMAT(length(SXMAT)+1) = SX;
        SYMAT(length(SYMAT)+1) = SY;
        SXYMAT(length(SXYMAT)+1) = SXY;
    end
end
end

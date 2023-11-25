function [SFL,GDSFL,JAC] = INTERPLN2D(NPE,XI,ETA,ELXY)

NP = [1 2 3 4 5 7 6 8 9];
XNODE = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0];

if NPE  == 4
    for I = 1:NPE
        XP = XNODE(I,1);
        YP = XNODE(I,2);
        XI0 = 1 + XI*XP;
        ETA0 = 1 + ETA*YP;

        SFL(I) = 0.25*XI0*ETA0;
        DSFL(1,I) = 0.25*XP*ETA0;
        DSFL(2,I) = 0.25*YP*XI0;
    end
elseif NPE == 8
    for I = 1:NPE
        NI = NP(I);
        XP = XNODE(NI,1);
        YP = XNODE(NI,2);
        XI0 = 1 + XP*XI;
        ETA0 = 1 + YP*ETA;
        XI1 = 1 - XI^2;
        ETA1 = 1 - ETA^2;
        if I<=4
            SFL(NI) = 0.5*XI0*ETA0*(XI*XP + ETA*YP - 1);
            DSFL(1,NI) = 0.25*ETA0*XP*(2*XI*XP + ETA*YP);
            DSFL(2,NI) = 0.25*XI0*YP*(2*ETA*YP + XI*XP);
        else
            if I <= 6
                SFL(NI) = 0.5*XI1*ETA0;
                DSFL(1,NI) = -1*XI*ETA0;
                DSFL(2,NI) = 0.5*YP*XI1;
            else
                SFL(NI) = 0.5*ETA1*XI0;
                DSFL(1,NI) = 0.5*XP*ETA1;
                DSFL(2,NI) = -ETA*XI0;
            end
        end
    end
elseif NPE == 9
    for I = 1:NPE
        NI = NP(I);
        XP = XNODE(NI,1);
        YP = XNODE(NI,2);
        XI0 = 1 + XP*XI;
        ETA0 = 1 + YP*ETA;
        XI1  = 1.0-XI*XI;
        ETA1 = 1.0-ETA*ETA;
        XI2 = XP*XI;
        ETA2 = YP*ETA;
        
        if I <= 4
            SFL(NI) = 0.25*XI0*ETA0*XI2*ETA2;
            DSFL(1,NI) = 0.25*XP*ETA2*ETA0*(1 + 2*XI2);
            DSFL(2,NI) = 0.25*YP*XI2*XI0*(1 + 2*ETA2);
        else
            if I <= 6
                SFL(NI) = 0.5*XI1*ETA0*ETA2;
                DSFL(1,NI) = -XI*ETA2*ETA0;
                DSFL(2,NI) = 0.5*XI1*YP*(1 + 2*ETA2);
            else
                if I <=8
                    SFL(NI)    = 0.5*ETA1*XI0*XI2;
                    DSFL(2,NI) = -ETA*XI2*XI0;
                    DSFL(1,NI) = 0.5*ETA1*XP*(1.0+2.0*XI2);
                else
                    SFL(NI) = XI1*ETA1;
                    DSFL(1,NI) = -2*XI*ETA1;
                    DSFL(2,NI) = -2*ETA*XI1;
                end
            end
        end
    end
end

% Compute Jacobian Matrix

GJ = zeros(2);
for I = 1:2
    for J = 1:2
        for K = 1:NPE
            GJ(I,J) = GJ(I,J) + DSFL(I,K)*ELXY(K,J);
        end
    end
end

GJINV = GJ^-1;  % Inverse of Jacobian Matrix

JAC = det(GJ);

% GDSFL  = GJINV * DSFL
GDSFL = zeros(2,NPE);
for I = 1:2
    for J = 1:NPE
        for K = 1:2
            GDSFL(I,J) = GDSFL(I,J) + GJINV(I,K)*DSFL(K,J);
        end
    end
end
end




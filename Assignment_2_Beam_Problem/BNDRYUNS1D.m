function [GLK,GLF] = BNDRYUNS1D(NONLIN,NDF,GLK,GLF,GLU,NSPV,ISPV,VSPV,NSSV,ISSV,VSSV,NSMB,ISMB,BETA0,BETAU,UREF)

%% Essential Boundary Condition
if NSPV > 0
    for NP = 1:NSPV
        NB = (ISPV(NP,1)-1)*NDF + ISPV(NP,2);
        if NONLIN == 1 %For Direct Iteration Method
            GLF(NB) = VSPV(NP);
        elseif NONLIN == 2 %For Newton's Method
            GLF(NB) = 0;
        end
        GLK(NB,:) = 0;
        GLK(NB,NB) = 1;
    end
end

%% Natural Boundary Condition
if NSSV>0
    for i = 1:NSSV
        NB = (ISSV(i,1)-1)*NDF+ISSV(i,2);
        GLF(NB)=GLF(NB)+VSSV(i);
    end
end

%%Mixed Boundary Condition
if(NSMB>0)
    for MB=1:NSMB
        NB = (ISMB(MB,1)-1)*NDF+ISMB(MB,2);
        if(NONLIN==1)
            GLK(NB,NB)=GLK(NB,NB)+BETA0(MB)+BETAU(MB)*GLU(NB);
            GLF(NB)=GLF(NB)+UREF(MB)*(BETA0(MB)+BETAU(MB)*GLU(NB));
        else
            GLK(NB,NB)= GLK(NB,NB) + BETA0(MB) + 2.0*BETAU(MB)*GLU(NB) - UREF(MB)*BETAU(MB);
            GLF(NB)= GLF(NB) + UREF(MB)*(BETA0(MB) + BETAU(MB)*GLU(NB)) - (BETA0(MB)+BETAU(MB)*GLU(NB))*GLU(NB);
        end
    end
end

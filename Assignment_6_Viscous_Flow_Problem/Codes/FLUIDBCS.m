function [GLK,GLF] = FLUIDBCS(NONLIN,NDF,NEQ,GLK,GLF,NSPV,ISPV,VSPV,NSSV,ISSV,VSSV)

%% Essential Boundary Condition

if NSPV > 0
    for NP = 1:NSPV
        NB = (ISPV(NP,1)-1)*NDF + ISPV(NP,2);
        for J = 1:NEQ
            GLK(NB,J) = 0;
            GLK(NB,NB) = 1;
            GLF(NB) = VSPV(NP);
        end
    end
end

%% Natural Boundary Condition

if NSSV > 0
    for NS = 1:NSSV
        NB = (ISSV(NS,1)-1)*NDF + ISSV(NS,2);
        GLF(NB) = GLF(NB) + VSSV(NS);
    end
end

end
function [GLK,GLF] = CONTBCS(NDF,NEQ,GLK,GLF,NSPV,ISPV,VSPV,NSSV,ISSV,VSSV)
%% Essential BCs
if NSPV~=0
    for NP = 1:NSPV
        NB = (ISPV(NP,1)-1)*NDF + ISPV(NP,2);

        GLK(NB,:) = 0;
        GLK(NB,NB) = 1;
        GLF(NB) = VSPV(NP);
    end
end

%% Natural BCs
if NSSV~=0
    for NS = 1:NSSV
        NB = (ISSV(NS,1)-1)*NDF + ISSV(NS,2);
        GLF(NB) = GLF(NB) + VSSV(NS);
    end
end

end

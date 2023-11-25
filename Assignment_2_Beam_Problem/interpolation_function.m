function [SFL,GDSFL,SFH,GDSFH,GDDSFH,GJ] = interpolation_function(ELX,IEL,model,NPE,xi)
% If a == 0, ie Lagrange interpolation function
% If a == 1, ie Hermit interpolation function
% Position of 1st element
% Position of last element
% xi is point of evaluation
syms x;
    
NPE = length(ELX);

he = ELX(end) - ELX(1);

if IEL == 1
    SFL = [0.5*(1-x); 0.5*(1+x)];
    DSFL = diff(SFL);
elseif IEL == 2
    NPE = IEL+1;
    SFL = [-0.5*x*(1-x); (1-x*x); 0.5*x*(1+x)];
    DSFL = diff(SFL);
end

% Jacobian
GJ = 0;
for i = 1:NPE
    GJ = GJ + DSFL(i)*ELX(i);
end

GDSFL = DSFL / GJ;

if model >=2
    SFH = [0.25*(2-3*x+x^3) ; -he*(1-x)*(1-x^2)/8 ; 0.25*(2+3*x-x^3) ; he*(1+x)*(1-x^2)/8];
    DSFH = diff(SFH);
    DDSFH = diff(DSFH);
    GDSFH = DSFH/GJ;
    GDDSFH = DDSFH/GJ/GJ;
end

x = xi;
SFL = eval(SFL);
GDSFL = eval(GDSFL);
SFH = eval(SFH);
GDSFH = eval(GDSFH);
GDDSFH = eval(GDDSFH);
end
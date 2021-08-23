function dydt = SEIR12EqnCalc(t,x,vp,fp)

% assigning values for variables from function inputs (vp,fp). 
N       = fp(1);
beta1   = fp(13);
beta2   = fp(14);
beta3   = fp(15);
beta4   = fp(16);
beta5   = fp(17);
lambda1 = fp(18);
lambda2 = fp(19);
lambda3 = fp(20);
theta1  = fp(21);
theta2  = fp(22);
tau     = fp(23);
tau2    = fp(24);
tau3    = fp(25);
tau4    = fp(26);


C1      = vp(1);
C2      = vp(2);
epsI    = vp(3);
eps0    = vp(4);
betaF   = vp(5);
muF     = vp(6);
f       = vp(7);
etaA    = vp(8);
etaH    = vp(9);
nu      = vp(10);
sigma   = vp(11);
phi     = vp(12);
gammaI  = vp(13);
gammaA  = vp(14);
gammaH  = vp(15);
gammaICU= vp(16);
deltaI  = vp(17);
deltaA  = vp(18);
deltaH  = vp(19);
deltaICU= vp(20);


Su      = x(1);
Sm      = x(2);
Eu      = x(3);
Em      = x(4);
Iu      = x(5);
Im      = x(6);
Au      = x(7);
Am      = x(8);
Ih      = x(9);
Iicu    = x(10);
R       = x(11);
D       = x(12);

% Implementing time dependent transmission coefficients. Refer to equation
% 16 and 17 in the manuscript. 
if t <= tau
    beta  = beta1+beta2;
elseif (t > tau)&&(t <= tau2)
    beta = beta1 + beta2*exp(-lambda1*(t-tau));
elseif (t > tau2)&&(t <= tau3)
     beta = beta1 + beta2*exp(-lambda1*(tau2-tau))+beta3*(1-exp(-lambda2*(t-tau2)));
elseif(t > tau3)&&(t <= tau4)
     beta = beta1 + beta2*exp(-lambda1*(tau2-tau))+beta3*(1-exp(-lambda2*(tau3-tau2)))+beta4*exp(-lambda1*(t-tau3)); 
else
     beta = beta1 + beta2*exp(-lambda1*(tau2-tau))+beta3*(1-exp(-lambda2*(tau3-tau2)))+beta4*exp(-lambda1*(tau4-tau3))+beta5*(1-exp(-lambda2*(t-tau4))); 
end

beta = beta*C1;


% implementing time dependent theta values. refer to equation 18 in the
% manuscript. 
if t < tau2
    theta = theta1;
else
    theta = theta2;
end
theta = theta*C2;

% implementation of equations 4-15 in the manuscript. 
Np = N-(D+Iicu);

dydt = zeros(12,1);
% dydt(1) = -             beta*Su*(Iu + etaA * Au + (1-eps0)*(Im + etaA * Am + etaH * Ih))/Np - betaF * (1 - exp(-lambda2 * D))*Su/Np + muF * Sm*(Su + R)/Np;
% dydt(2) = -(1 - epsI) * beta*Sm*(Iu + etaA * Au + (1-eps0)*(Im + etaA * Am + etaH * Ih))/Np + betaF * (1 - exp(-lambda2 * D))*Su/Np - muF * Sm*(Su + R)/Np;

dydt(1) = -             beta*Su*(Iu + etaA * Au + (1-eps0)*(Im + etaA * Am + etaH * Ih))/Np - betaF * (1 - exp(-lambda3 * D*t))*D*Su/Np + muF * Sm*(Su + R)/Np;
dydt(2) = -(1 - epsI) * beta*Sm*(Iu + etaA * Au + (1-eps0)*(Im + etaA * Am + etaH * Ih))/Np + betaF * (1 - exp(-lambda3 * D*t))*D*Su/Np - muF * Sm*(Su + R)/Np;
dydt(3) =               beta*Su*(Iu + etaA * Au + (1-eps0)*(Im + etaA * Am + etaH * Ih))/Np - sigma * Eu;
dydt(4) =  (1 - epsI) * beta*Sm*(Iu + etaA * Au + (1-eps0)*(Im + etaA * Am + etaH * Ih))/Np - sigma * Em;
dydt(5) =  (1 - f) * sigma * Eu - theta * Iu - phi * Iu - gammaI * Iu - deltaI * Iu;
dydt(6) =  (1 - f) * sigma * Em + theta * Iu - phi * Im - gammaI * Im - deltaI * Im;
dydt(7) =  f * sigma * Eu - gammaA * Au - deltaA * Au;
dydt(8) =  f * sigma * Em - gammaA * Am - deltaA * Am;
dydt(9) =  phi * Iu + phi * Im - gammaH * Ih - nu * Ih - deltaH * Ih;
dydt(10) = nu * Ih - gammaICU * Iicu - deltaICU * Iicu;
dydt(11) = gammaA * (Au + Am) + gammaI * (Iu + Im) + gammaH * Ih + gammaICU * Iicu;
dydt(12) = deltaA * (Au + Am) + deltaI * (Iu + Im) + deltaH * Ih + deltaICU * Iicu;


end

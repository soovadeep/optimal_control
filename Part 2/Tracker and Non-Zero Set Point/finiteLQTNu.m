function [nudot] = finiteLQTNu(t,nu)

global A B Rinv N Q SMat Amp w L 

zrdot = Amp*w*cos(w*t);
ref = [0.05; 0; 0; 0];

K = Rinv*B'*SMat;

nudot = SMat*B*Rinv*B'*nu + SMat*B*Rinv*N'*ref + SMat*L*zrdot - A'*nu - Q*ref + N*Rinv*B'*nu + N*Rinv*N'*ref; 
end

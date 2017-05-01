function [y,zsddot,F] = car_lqt_finite(t,x)

global ks bs ms w Amp A B L SMat Rinv N nuiter 

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

ref = [0.05; 0; 0; 0];

K3 = Rinv*(B'*SMat + N');
F = -K3*[x(1);x(2);x(3);x(4)] + Rinv*N'*ref + Rinv*B'*nuiter;

zrdot = Amp*w*cos(w*t);

y = A*x + B*F + L*zrdot;
zsddot = (F - ks*x1 - bs*(x2 - x4))/ms;
end

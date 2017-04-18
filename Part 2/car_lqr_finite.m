function [y,zsddot] = car_lqr_finite(t,x)

global ks bs ms w Amp A B L SMat Rinv N

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

K2 = Rinv*(B'*SMat + N');
F = -K2*x;

zrdot = Amp*w*cos(w*t);

y = A*x + B*F + L*zrdot;
zsddot = (F - ks*x1 - bs*(x2 - x4))/ms;
end

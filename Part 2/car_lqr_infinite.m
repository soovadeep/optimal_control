function [y,zsddot] = car_lqr_infinite(t,x)

global ks bs ms w Amp A B L K

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

F = -K*x;

zrdot = Amp*w*cos(w*t);

y = A*x + B*F + L*zrdot;
zsddot = (F - ks*x1 - bs*(x2 - x4))/ms;
end

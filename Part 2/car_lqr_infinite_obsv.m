function [xdot,zsddot] = car_lqr_infinite_obsv(t,x)

global ks bs ms w Amp A B L K FW G H C1 M

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
W = x(5:6);
ytilde = x(7:8);

z = C1*x(1:2);
yhat = M*z + W;

F = -K*[x1;x2;yhat(1);yhat(2)];

zrdot = Amp*w*cos(w*t);

statedot = A*x(1:4) + B*F + L*zrdot;
wdot = FW*W + G*z + H*F;
ytildedot = FW*ytilde;
xdot = [statedot;wdot;ytildedot];
zsddot = (F - ks*x1 - bs*(x2 - x4))/ms;
end

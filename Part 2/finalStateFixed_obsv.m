function [xdot,zsddot, F] = finalStateFixed_obsv(t,x)

global A B Rinv SMat VVec PVec ks bs ms w Amp N FW G H C1 M L

% SMat = reshape(S, [4,4]);
% VVec = reshape(V, [4,4]);
% PVec = reshape(P, [4,4]);

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
W = x(5:6);
ytilde = x(7:8);

z = C1*x(1:2);
yhat = M*z + W;

rFinal = zeros(1,1);
K = Rinv*(B'*SMat + N');

F = -(K - Rinv*B'*VVec*inv(PVec)*VVec')*[x1;x2;yhat(1);yhat(2)] - Rinv*B'*VVec*inv(PVec)*rFinal; % Control law


% xdot = A*x + B*u;

zrdot = Amp*w*cos(w*t);
statedot = A*x(1:4) + B*F + L*zrdot;
wdot = FW*W + G*z + H*F;
ytildedot = FW*ytilde;
xdot = [statedot;wdot;ytildedot];
zsddot = (F - ks*x1 - bs*(x2 - x4))/ms;

end
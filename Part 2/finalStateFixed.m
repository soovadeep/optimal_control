function [xdot,zsddot, u] = finalStateFixed(t,x)

global A B Rinv SMat VVec PVec ks bs ms w Amp N

% SMat = reshape(S, [4,4]);
% VVec = reshape(V, [4,4]);
% PVec = reshape(P, [4,4]);

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

rFinal = 0.05;
K = Rinv*(B'*SMat + N');

u = -(K - Rinv*B'*VVec*inv(PVec)*VVec')*x - Rinv*B'*VVec*inv(PVec)*rFinal;

xdot = A*x + B*u;
zsddot = (u - ks*x1 - bs*(x2 - x4))/ms;

end
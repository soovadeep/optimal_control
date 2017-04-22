function [xdot] = finalStateFixed(x, S, V, P)

global A B Rinv

SMat = reshape(S, [4,4]);
VMat = reshape(V, [4,4]);
PMat = reshape(P, [4,4]);

rFinal = zeros(4,1);
K = Rinv*B*SMat;

u = -(K - Rinv*B'*VMat*inv(PMat)*VMat')*x - Rinv*B'*VMat*inv(PMat)*rFinal;

xdot = A*x + B*u;

end
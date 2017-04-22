function [Vdot] = finalStateFixed_V(S,V)
global A B Rinv 

SMat = reshape(S, [4,4]);
VMat = reshape(V, [4,4]);

VdotMat = (-A' + SMat*B*Rinv*B'*SMat)*VMat;

Vdot = reshape(VdotMat', [16,1]);
end
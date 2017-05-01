function [Vdot] = finalStateFixed_V(t,V)
global Abar B Rinv SMat

% SMat = reshape(S, [4,4]);
% VMat = reshape(V, [4,4]);

Vdot = (-Abar' + SMat*B*Rinv*B')*V;

% Vdot = reshape(VdotMat', [16,1]);
end
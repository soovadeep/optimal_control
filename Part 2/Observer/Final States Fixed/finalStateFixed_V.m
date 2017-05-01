function [Vdot] = finalStateFixed_V(t,V)
global Abar B Rinv SMat

Vdot = (-Abar' + SMat*B*Rinv*B')*V;

end
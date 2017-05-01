function [Pdot] = finalStateFixed_P(t,P)

global B Rinv VVec

Pdot = (VVec'*B*Rinv*B'*VVec);
end

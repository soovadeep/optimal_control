function [Pdot] = finalStateFixed_P(t,P)

global B Rinv VVec

% VMat = reshape(V,[4,4]);

Pdot = (VVec'*B*Rinv*B'*VVec);

% Pdot = reshape(PdotMat', [16, 1]);

end

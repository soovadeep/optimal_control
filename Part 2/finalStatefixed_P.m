function [Pdot] = finalStatefixed_P(V,P)

global B Rinv

VMat = reshape(V,[4,4]);

PdotMat = (VMat'*B*Rinv*B'*VMat);

Pdot = reshape(PdotMat', [16, 1]);

end

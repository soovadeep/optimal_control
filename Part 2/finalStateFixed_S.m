function [Sdot] = finalStateFixed_S(S)
global A B Rinv Q

SdotMat = -(S*A + A'*S - S*B*Rinv*B'*S + Q);
Sdot =  reshape(SdotMat',[16,1]);
end
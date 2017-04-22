function [Sdot] = finalStateFixed_S(t,S)
global A B Rinv Q Qbar Abar N

SMat = reshape(S,[4,4])';

% SdotMat = -(S*A + A'*S - S*B*Rinv*B'*S + Q);
% SdotMat = -Abar'*SMat - SMat*Abar + SMat*B*Rinv*B'*SMat - Qbar;
SdotMat = -((A - B*Rinv*N')'*SMat + SMat*(A - B*Rinv*N') + (Q - N*Rinv*N') - SMat*B*Rinv*B'*SMat);
Sdot =  reshape(SdotMat',[16,1]);
end
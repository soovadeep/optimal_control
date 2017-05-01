function [Sdot] = finiteLQTRiccati(t,S)

global A B Rinv N Q Qbar Abar

SMat = reshape(S,[4,4])';

SdotMat = -((A - B*Rinv*N')'*SMat + SMat*(A - B*Rinv*N') + (Q - N*Rinv*N') - SMat*B*Rinv*B'*SMat);
% SdotMat = -Abar*SMat - SMat*Abar + SMat*B*Rinv*B'*Smat - Qbar;
Sdot =  reshape(SdotMat',[16,1]);
end

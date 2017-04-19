function [Sdot] = finiteLQTRiccati(t,S)

global A B Rinv N Q

SMat = reshape(S,[4,4])';

SdotMat = -((A - B*Rinv*N')'*SMat + SMat*(A - B*Rinv*N') + (Q - N*Rinv*N') - SMat*B*Rinv*B'*SMat);
Sdot =  reshape(SdotMat',[16,1]);
end

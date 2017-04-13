function residual = OLcontrolBC(ya,yb)

global x10 tf Amp w A L Q

zr = Amp*sin(w*tf); 
zrdot = Amp*w*cos(w*tf); 

% X = [-zr;0;0;zrdot];
% Lambda = -1/2*((A*X + L*zrdot)^-1)'*(X'*Q*X);

residual = [(ya(1) - x10); ya(2); ya(3); ya(4); (yb(5)); (yb(6)); (yb(7)); (yb(8))];
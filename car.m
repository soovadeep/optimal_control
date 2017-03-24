function y=car(t,x)

global w Amp A B L 

u = 0;
zrdot = Amp*w*cos(w*t);

y = A*x+B*u+L*zrdot;

end

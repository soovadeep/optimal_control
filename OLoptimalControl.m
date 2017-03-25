function [xdot, F, zsddot]=OLoptimalControl(t,x) 

global ks ms bs w Amp A B L N Q Rinv

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

X = [x1; x2; x3; x4];

zrdot = Amp*w*cos(w*t);

lambda1 = x(5);
lambda2 = x(6);
lambda3 = x(7);
lambda4 = x(8);

lambda = [lambda1; lambda2; lambda3; lambda4];

% Stationarity equation

F = -Rinv*B'*lambda;

% State and Co-state equations

Xdot = A*X - (B*Rinv*B')*lambda + L*zrdot -(B*Rinv)*X'*N;
lambdadot = -Q*X + (N*Rinv*B' - A')*lambda + N*Rinv*X'*N;

xdot = [Xdot;lambdadot];

zsddot = (F - ks*x1 - bs*(x2 - x4))/ms;
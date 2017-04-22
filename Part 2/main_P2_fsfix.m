clear

global ks kt ms mu bs bt w Amp A B C L rho1 rho2 rho3 rho4 R N Q Rinv x10 x20 x30 x40 tf K SMat t0 nuiter Qbar Abar

kt = (704*10^3)/4; % N/m
ks = 15*10^3; % N/m
bs = 1400; % Ns/m  Check
bt = 0; % Ns/m
mu = 181/4; % kg
ms = 1814/4; % kg
rho1 = 0.4; % 0.4 
rho2 = 0.04; % 0.04
rho3 = 0.4; % 0.4
rho4 = 0.04; % 0.04
Amp = 0.05;
w = 0*2*pi;
t0 = 0;
tf = 10;
steps = tf*1000;
stepsize = (tf-t0)/steps;

A = [0 1 0 -1; -ks/ms -bs/ms 0 bs/ms; 0 0 0 1; ks/mu bs/mu -kt/mu -(bs+bt)/mu ];
B = [0;1/ms;0;-1/mu];
L = [0;0;-1;0];
C = [1 0 0 0;0 1 0 0];

R = 1/ms^2;
Rinv = 1/R;
N = [-ks/ms^2; -bs/ms^2; 0; bs/ms^2];
Q = [(ks^2/ms^2 + rho1)  bs*ks/ms^2            0      -bs*ks/ms^2;
     bs*ks/ms^2          (bs^2/ms^2 + rho2)    0      -bs^2/ms^2;
     0                   0                     rho3   0;
     -bs*ks/ms^2         -bs^2/ms^2            0      (bs^2/ms^2 + rho2)];

zs0 = -0.05;
zu0 = 0;
zsdot0 = 0;
zudot0 = 0;
zr0 = 0;

x10 = zs0-zu0;
x20 = zsdot0;
x30 = zu0-zr0;
x40 = zudot0;

x0 = [x10; x20; x30; x40];

tspan = [t0 tf];
tspan_ulta = [tf t0];
ST = zeros(4,4);

[T, S] = rk4fixed(@finalStateFixed_S, tspan_ulta, ST, steps);

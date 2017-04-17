close all
clear
clc

%% Parameters and Passive Data

global ks kt ms mu bs bt w Amp A B C L rho1 rho2 rho3 rho4 R N Q Rinv x10 x20 x30 x40 tf K S

load('Y_passive_dynamic.mat')
Tpass = T;
Ypass = Y;

clear T Y;

kt = (704*10^3)/4; % N/m
ks = 15*10^3; % N/m
bs = 1400; % Ns/m  Check
bt = 0; % Ns/m
mu = 181/4; % kg
ms = 1814/4; % kg
rho1 = 0.4; 
rho2 = 0.04;
rho3 = 0.4;
rho4 = 0.04;
Amp = 0.05;
w = 2*pi;
t0 = 0;
tf = 15;
steps = 15000;

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

%% Infinite time LQR with no road input

% x1 = zs-zu
% x2 = zsdot
% x3 = zu-zr
% x4 = zudot

[K,S,e] = lqr(A,B,Q,R,N);

zs0 = 0;
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
[T,Y] = rk4fixed(@car_lqr_infinite,tspan,x0,steps);

lengthPass = size(T,1);
zacclPass = zeros(1,lengthPass);

for i = 1:lengthPass
    [xdot, zsddotPass] = car(T(i),Y(i,:)');
    zacclPass(1,i) = zsddotPass;
end

u = 0;
ZR = Amp*sin(w*T);

zu = Y(:,3) + ZR;
zs = Y(:,1) + zu;

%%
fig = figure(1);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,zs,'-r','LineWidth',1.5)
hold on
plot(T,ZR,'-.b','LineWidth',1.5)
title('Sprung Mass Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
legend('Response','Road Profile')
set(legend,'Interpreter','Latex','FontSize',12)
% print('Passive-SMD','-djpeg','-r300')

fig = figure(2);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,zacclPass,'-r','LineWidth',1.5)
title('Sprung Mass Acceleration vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)
% print('Passive-SMA','-djpeg','-r300')

fig = figure(3);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,Y(:,1),'-r','LineWidth',1.5)
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-SD','-djpeg','-r300')

fig = figure(4);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,Y(:,3),'-r','LineWidth',1.5)
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-TD','-djpeg','-r300')
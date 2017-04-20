% close all
clear
% clc

%% Parameters and Passive Data

global ks kt ms mu bs bt w Amp A B C L rho1 rho2 rho3 rho4 R N Q Rinv x10 x20 x30 x40 tf K SMat t0 nuiter Qbar Abar FW H G M C1

load('Y_passive_static.mat')
TPass = T;
YPass = Y;
zuPass = zu;
zsPass = zs;

clear T Y zu zs;

load('Y_infinitelqr_static.mat')

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
A11 = A(1:2,1:2);
A12 = A(1:2,3:4);
A21 = A(3:4,1:2);
A22 = A(3:4,3:4);
B = [0;1/ms;0;-1/mu];
B1 = B(1:2);
B2 = B(3:4);
L = [0;0;-1;0];
L2 = L(3:4);
C = [1 0 0 0;0 1 0 0];
C1 = C(1:2,1:2);

R = 1/ms^2;
Rinv = 1/R;
N = [-ks/ms^2; -bs/ms^2; 0; bs/ms^2];
Q = [(ks^2/ms^2 + rho1)  bs*ks/ms^2            0      -bs*ks/ms^2;
     bs*ks/ms^2          (bs^2/ms^2 + rho2)    0      -bs^2/ms^2;
     0                   0                     rho3   0;
     -bs*ks/ms^2         -bs^2/ms^2            0      (bs^2/ms^2 + rho2)];
 

[K,~,e] = lqr(A,B,Q,R,N);

p = [real(25*e(1));real(20*e(1))];
M = place(A22',(C1*A12)',p);

FW = A22 - M*C1*A12;
H = B2 - M*C1*B1;
G = (A21 - M*C1*A11)*inv(C1) + FW*M;

zs0 = -0.05;
zu0 = 0;
zsdot0 = 0;
zudot0 = 0;
zr0 = 0;

x10 = zs0-zu0;
x20 = zsdot0;
x30 = zu0-zr0;
x40 = zudot0;
x50 = -0.01;
x60 = -0.01;
x70 = -0.01;
x80 = -0.01;

x0 = [x10; x20; x30; x40; x50; x60; x70; x80];

tspan = [t0 tf];
[Tobsv,Yobsv] = rk4fixed(@car_lqr_infinite_obsv,tspan,x0,steps);

lengthPass = size(Tobsv,1);
zacclObsv = zeros(1,lengthPass);

for i = 1:lengthPass
    [xdot, zsddot] = car_lqr_infinite_obsv(Tobsv(i),Yobsv(i,:)');
    zacclObsv(1,i) = zsddot;
end

u = 0;
ZR = Amp*sin(w*Tobsv);

zuObsv = Yobsv(:,3) + ZR;
zsObsv = Yobsv(:,1) + zuObsv;

fig = figure(10);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(Tobsv,zsObsv,'-g','LineWidth',1.5)
hold on
plot(T(1:steps),zs(1:steps),'-k','LineWidth',1.5)
plot(TPass(1:steps),zsPass(1:steps),'-r','LineWidth',1.5)
plot(Tobsv,ZR,'-.b','LineWidth',1.5)
title('Sprung Mass Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
legend('Active (Observer)','Active','Passive','Road Profile')
set(legend,'Interpreter','Latex','FontSize',12)
% print('Passive-SMD','-djpeg','-r300')

fig = figure(11);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(Tobsv,zacclObsv,'-g','LineWidth',1.5)
hold on
plot(T(1:steps),zaccl(1:steps),'-k','LineWidth',1.5)
plot(TPass(1:steps),zacclPass(1:steps),'-r','LineWidth',1.5)
legend('Active (Observer)','Active','Passive')
title('Sprung Mass Acceleration vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)
% print('Passive-SMA','-djpeg','-r300')

fig = figure(12);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(Tobsv,Yobsv(:,1),'-g','LineWidth',1.5)
hold on
plot(T(1:steps),Y(1:steps,1),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,1),'-r','LineWidth',1.5)
legend('Active (Observer)','Active','Passive')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-SD','-djpeg','-r300')

fig = figure(13);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(Tobsv,Yobsv(:,3),'-g','LineWidth',1.5)
hold on
plot(T(1:steps),Y(1:steps,3),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
legend('Active (Observer)','Active','Passive')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-TD','-djpeg','-r300')




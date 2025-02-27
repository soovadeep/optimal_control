close all
clear
% clc

%% Parameters and Passive Data

global ks kt ms mu bs bt w Amp A B C L rho1 rho2 rho3 rho4 R N Q Rinv x10 x20 x30 x40 tf K SMat t0 FW H G M C1 VVec PVec

load('Y_passive_static.mat')
TPass = T;
YPass = Y;
zuPass = zu;
zsPass = zs;

clear T Y zu zs;

load('Y_finitelqr_fsf_static.mat')
YFSF = Y;
TFSF = T; 
zuFSF = zu;
zsFSF = zs;
zacclFSF = zaccl;
forceFSF = force;

clear T Y zu zs zaccl force;

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

%% Infinite LQR w/ Observer

p = [real(20*e(1));real(25*e(1))];
M = place(A22',(C1*A12)',p);

FW = A22 - M*C1*A12;
H = B2 - M*C1*B1;
G = (A21 - M*C1*A11)*inv(C1) + FW*M;

%{
zs0 = 0;
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
forceObsv = zeros(1,lengthPass);

for i = 1:lengthPass
    [xdot, zsddot,F] = car_lqr_infinite_obsv(Tobsv(i),Yobsv(i,:)');
    zacclObsv(1,i) = zsddot;
    forceObsv(1,i) = F;
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
% set(legend,'Interpreter','Latex','FontSize',12)
print('Passive-SMD-Observer-Infinite-Dynamic','-djpeg','-r300')

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
print('Passive-SMA-Observer-Infinite-Dynamic','-djpeg','-r300')

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
print('Passive-SD-Observer-Infinite-Dynamic','-djpeg','-r300')

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
print('Passive-TD-Observer-Infinite-Dynamic','-djpeg','-r300')

fig = figure(20);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(Tobsv,forceObsv,'-g','LineWidth',1.5)
hold on 
plot(T,force,'-k','LineWidth',1.5)
% plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
legend('Active (Observer)','Active')
title('Control Input vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Force\hspace{0.05in}(N)$','Interpreter','Latex','FontSize',12)
print('Passive-Force-Observer-Infinite-Dynamic','-djpeg','-r300')
%}

%% Finite LQR w/ Observer

%{
clear S;

tspan = [tf t0-stepsize];
S0 = zeros(16,1);
[~,S] = rk4fixed(@finiteLQRRiccati,tspan,S0,steps+1);

S = flipud(S);

zs0FTObsv = -0.05;
zu0FTObsv = 0;
zsdot0FTObsv = 0;
zudot0FTObsv = 0;
zr0FTObsv = 0;

x10FTObsv = zs0FTObsv-zu0FTObsv;
x20FTObsv = zsdot0FTObsv;
x30FTObsv = zu0FTObsv-zr0FTObsv;
x40FTObsv = zudot0FTObsv;
x50FTObsv = -0.01;
x60FTObsv = -0.01;
x70FTObsv = -0.01;
x80FTObsv = -0.01;

x0FTObsv = [x10FTObsv; x20FTObsv; x30FTObsv; x40FTObsv; x50FTObsv; x60FTObsv; x70FTObsv; x80FTObsv];
YFTObsv = zeros(steps,8);
TFTObsv = zeros(steps,1);
zacclFTObsv = zeros(steps,1);
forceFTObsv = zeros(steps,1);
YFTObsv(1,:) = x0FTObsv';
K2Obsv = zeros(steps,4);

t0iter = t0;

for i = 1:steps-1
    i
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    tfiter = t0iter + stepsize;
    tspan = [t0iter tfiter];
    x0FTiter = YFTObsv(i,:)';
    [TFTiterObsv,YFTiterObsv] = ode15s(@car_lqr_finite_obsv,tspan,x0FTiter); %400
    t0iter = tfiter;
    K2Obsv(i,:) = Rinv*(B'*SMat + N');
    [YdotFTObsv, zsddotFTObsv, FFT] = car_lqr_finite_obsv(TFTObsv(i),YFTObsv(i,:)');
%     YFT(i+1,:) = YFT(i,:) + YdotFT'*stepsize;
    zacclFTObsv(i) = zsddotFTObsv;
    forceFTObsv(i) = FFT;
    TFTObsv(i+1) = tfiter;
    YFTObsv(i+1,:) = YFTiterObsv(end,:);
end

ZR = Amp*sin(w*TFTObsv);

zuFTObsv = YFTObsv(:,3) + ZR;
zsFTObsv = YFTObsv(:,1) + zuFTObsv;

%%
fig = figure(18);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
% plot(Tobsv,zsObsv,'-g','LineWidth',1.5)
hold on
plot(TFTObsv,zsFTObsv,'-g','LineWidth',1.5)
plot(TFT(1:steps),zsFT(1:steps),'-k','LineWidth',1.5)
plot(TPass(1:steps),zsPass(1:steps),'-r','LineWidth',1.5)
plot(TFTObsv,ZR,'-.b','LineWidth',1.5)
title('Sprung Mass Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% legend('Active (Observer)','Active (Observer, Finite)','Active','Passive','Road Profile')
legend('Active (Observer, Finite)','Active','Passive','Road Profile')
% set(legend,'Interpreter','Latex','FontSize',12)
print('Passive-SMD-FLQR-Observer','-djpeg','-r300')

fig = figure(19);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
% plot(Tobsv,zacclObsv,'-g','LineWidth',1.5)
hold on
plot(TFTObsv,zacclFTObsv,'-g','LineWidth',1.5)
plot(TFT(1:steps),zacclFT(1:steps),'-k','LineWidth',1.5)
plot(TPass(1:steps),zacclPass(1:steps),'-r','LineWidth',1.5)
% legend('Active (Observer)','Active (Observer, Finite)','Active','Passive')
legend('Active (Observer, Finite)','Active','Passive')
title('Sprung Mass Acceleration vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)
print('Passive-SMA-FLQR-Observer','-djpeg','-r300')

fig = figure(20);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
% plot(Tobsv,Yobsv(:,1),'-g','LineWidth',1.5)
hold on
plot(TFTObsv,YFTObsv(:,1),'-m','LineWidth',1.5)
plot(TFT(1:steps),YFT(1:steps,1),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,1),'-r','LineWidth',1.5)
% legend('Active (Observer)','Active (Observer, Finite)','Active','Passive')
legend('Active (Observer, Finite)','Active','Passive')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
print('Passive-SD-FLQR-Observer','-djpeg','-r300')

fig = figure(21);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
% plot(Tobsv,Yobsv(:,3),'-g','LineWidth',1.5)
hold on
plot(TFTObsv,YFTObsv(:,3),'-m','LineWidth',1.5)
plot(TFT(1:steps),YFT(1:steps,3),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
% legend('Active (Observer)','Active (Observer, Finite)','Active','Passive')
legend('Active (Observer, Finite)','Active','Passive')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
print('Passive-TD-FLQR-Observer','-djpeg','-r300')

fig = figure(22);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TFTObsv,forceFTObsv,'-g','LineWidth',1.5)
hold on 
plot(TFT,forceFT,'-k','LineWidth',1.5)
% plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
legend('Active (Observer)','Active')
title('Control Input vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Force\hspace{0.05in}(N)$','Interpreter','Latex','FontSize',12)
print('Passive-Force-FLQR-Observer','-djpeg','-r300')
%}

%% LQR w/ FSF w/ Observer

tspan = [tf t0-stepsize];
S0 = zeros(16,1);
[~, S] = rk4fixed(@finalStateFixed_S, tspan, S0, steps+1);

S = flipud(S);

%%
V0 = [1 0 0 0]';
V = zeros(steps+1,4);
V(end,:) = V0';

tfiter = tf;

for i = steps + 1:-1:2
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    t0iter = tfiter - stepsize;
    tspan = [tfiter t0iter];
    V0iter = V(i,:)';
    [~,Viter] = ode23s(@finalStateFixed_V,tspan,V0iter);
    tfiter = t0iter;
    V(i-1,:) = Viter(end,:);
end

P0 = [0];
P = zeros(steps+1,1);
P(end,:) = P0';

tfiter = tf;

for i = steps + 1:-1:2
    VVec = V(i,:)';
    t0iter = tfiter - stepsize;
    tspan = [tfiter t0iter];
    P0iter = P(i)';
    [~,Piter] = ode15s(@finalStateFixed_P,tspan,P0iter);
    tfiter = t0iter;
    P(i-1) = Piter(end);
end

zs0 = -0.05;
zu0 = 0;
zsdot0 = 0;
zudot0 = 0;
zr0 = 0;

x10 = zs0-zu0;
x20 = zsdot0;
x30 = zu0-zr0;
x40 = zudot0;
x50FTObsv = -0.01;
x60FTObsv = -0.01;
x70FTObsv = -0.01;
x80FTObsv = -0.01;

x0 = [x10; x20; x30; x40; x50FTObsv; x60FTObsv; x70FTObsv; x80FTObsv];

Y = zeros(steps,8);
T = zeros(steps,1);
zaccl = zeros(steps,1);
force = zeros(steps,1);
Y(1,:) = x0';
K2 = zeros(steps,4);

t0iter = t0;

for i = 1:steps-1
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    VVec = V(i,:)';
    PVec = P(i);
    tfiter = t0iter + stepsize;
    tspan = [t0iter tfiter];
    x0iter = Y(i,:)';
    [Titer,Yiter] = ode15s(@finalStateFixed_obsv,tspan,x0iter); %400
    t0iter = tfiter;
    K2(i,:) = Rinv*(B'*SMat + N');
    [Ydot, zsddot, u] = finalStateFixed_obsv(T(i),Y(i,:)');
    zaccl(i) = zsddot;
    force(i) = u;
    T(i+1) = tfiter;
    Y(i+1,:) = Yiter(end,:);
end

ZR = Amp*sin(w*T);

zu = Y(:,3) + ZR;
zs = Y(:,1) + zu;

steps = tf*1000;

fig = figure(1);
clear title
clear legend
plot(T,Y(:,1),'-g','LineWidth',1.5)
hold on
plot(TPass(1:steps),YPass(1:steps),'-r','LineWidth',1.5)
plot(TFSF(1:steps),YFSF(1:steps),'-.k','LineWidth',1.5)
plot(T,ZR,'-.b','LineWidth',1.5)
title('Sprung Mass Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
legend('Active (FSF, Observer)','Passive','Active (FSF)','Road Profile')

fig = figure(2);
clear title
clear legend
plot(T,zaccl,'-g','LineWidth',1.5)
hold on
plot(TPass(1:steps),zacclPass(1:steps),'-r','LineWidth',1.5)
plot(TFSF(1:steps),zacclFSF(1:steps),'-.k','LineWidth',1.5)
legend('Active (FSF, Observer)','Passive','Active (FSF)')
title('Sprung Mass Acceleration vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)

fig = figure(20);
clear title
clear legend
hold on
plot(T,Y(:,1),'-m','LineWidth',1.5)
plot(TFSF(1:steps),YFSF(1:steps,1),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,1),'-r','LineWidth',1.5)
legend('Active (FSF, Observer)','Active (FSF)','Passive')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)

fig = figure(21);
clear title
clear legend
hold on
plot(T,Y(:,3),'-m','LineWidth',1.5)
plot(TFSF(1:steps),YFSF(1:steps,3),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
legend('Active (FSF, Observer)','Active (FSF)','Passive')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)

fig = figure(22);
clear title
clear legend
plot(T,force,'-g','LineWidth',1.5)
hold on 
plot(TFSF,forceFSF,'-k','LineWidth',1.5)
legend('Active (FSF, Observer)','Active (FSF)')
title('Control Input vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Force\hspace{0.05in}(N)$','Interpreter','Latex','FontSize',12)
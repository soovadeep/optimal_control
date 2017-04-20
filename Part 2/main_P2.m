close all
clear
clc

%% Parameters and Passive Data

global ks kt ms mu bs bt w Amp A B C L rho1 rho2 rho3 rho4 R N Q Rinv x10 x20 x30 x40 tf K SMat t0 nuiter Qbar Abar

load('Y_passive_dynamic.mat')
TPass = T;
YPass = Y;
zuPass = zu;
zsPass = zs;

clear T Y zu zs;

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
tf = 5;
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

%% Infinite time LQR

%{

% x1 = zs-zu
% x2 = zsdot
% x3 = zu-zr
% x4 = zudot

[K,~,e] = lqr(A,B,Q,R,N);

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
[T,Y] = rk4fixed(@car_lqr_infinite,tspan,x0,steps);

lengthPass = size(T,1);
zaccl = zeros(1,lengthPass);

for i = 1:lengthPass
    [xdot, zsddot] = car_lqr_infinite(T(i),Y(i,:)');
    zaccl(1,i) = zsddot;
end

u = 0;
ZR = Amp*sin(w*T);

zu = Y(:,3) + ZR;
zs = Y(:,1) + zu;

%}

%{
fig = figure(1);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,zs,'-g','LineWidth',1.5)
hold on
plot(TPass,zsPass,'-r','LineWidth',1.5)
plot(T,ZR,'-.b','LineWidth',1.5)
title('Sprung Mass Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
legend('Active','Passive','Road Profile')
set(legend,'Interpreter','Latex','FontSize',12)
% print('Passive-SMD','-djpeg','-r300')

fig = figure(2);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,zaccl,'-g','LineWidth',1.5)
hold on
plot(TPass,zacclPass,'-r','LineWidth',1.5)
legend('Active','Passive')
title('Sprung Mass Acceleration vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)
% print('Passive-SMA','-djpeg','-r300')

fig = figure(3);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,Y(:,1),'-g','LineWidth',1.5)
hold on
plot(TPass,YPass(:,1),'-r','LineWidth',1.5)
legend('Active','Passive')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-SD','-djpeg','-r300')

fig = figure(4);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,Y(:,3),'-g','LineWidth',1.5)
hold on 
plot(TPass,YPass(:,3),'-r','LineWidth',1.5)
legend('Active','Passive')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-TD','-djpeg','-r300')
%}

%% Finite time LQR

%{
clear S;

tspan = [tf t0-stepsize];
S0 = zeros(16,1);
[~,S] = rk4fixed(@finiteLQRRiccati,tspan,S0,steps+1);

S = flipud(S);

zs0FT = -0.05;
zu0FT = 0;
zsdot0FT = 0;
zudot0FT = 0;
zr0FT = 0;

x10FT = zs0FT-zu0FT;
x20FT = zsdot0FT;
x30FT = zu0FT-zr0FT;
x40FT = zudot0FT;

x0FT = [x10FT; x20FT; x30FT; x40FT];
YFT = zeros(steps,4);
TFT = zeros(steps,1);
zacclFT = zeros(steps,1);
YFT(1,:) = x0FT';
K2 = zeros(steps,4);

t0iter = t0;

for i = 1:steps-1
    i
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    tfiter = t0iter + stepsize;
    tspan = [t0iter tfiter];
    x0FTiter = YFT(i,:)';
    [TFTiter,YFTiter] = ode15s(@car_lqr_finite,tspan,x0FTiter); %400
    t0iter = tfiter;
    K2(i,:) = Rinv*(B'*SMat + N');
    [YdotFT, zsddotFT] = car_lqr_finite(TFT(i),YFT(i,:)');
%     YFT(i+1,:) = YFT(i,:) + YdotFT'*stepsize;
    zacclFT(i) = zsddotFT;
    TFT(i+1) = tfiter;
    YFT(i+1,:) = YFTiter(end,:);
end

zuFT = YFT(:,3) + ZR;
zsFT = YFT(:,1) + zuFT;

fig = figure(1);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TFT,zsFT,'-g','LineWidth',1.5)
hold on
plot(T,zs,'-k','LineWidth',1.5)
plot(TPass(1:steps),zsPass(1:steps),'-r','LineWidth',1.5)
plot(T,ZR,'-.b','LineWidth',1.5)
title('Sprung Mass Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
legend('Active (FT)','Active (IT)','Passive','Road Profile')
set(legend,'Interpreter','Latex','FontSize',12)
% print('Passive-SMD','-djpeg','-r300')

fig = figure(2);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TFT,zacclFT,'-g','LineWidth',1.5)
hold on
plot(T,zaccl,'-k','LineWidth',1.5)
plot(TPass(1:steps),zacclPass(1:steps),'-r','LineWidth',1.5)
legend('Active (FT)','Active (IT)','Passive')
title('Sprung Mass Acceleration vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)
% print('Passive-SMA','-djpeg','-r300')

fig = figure(3);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TFT,YFT(:,1),'-g','LineWidth',1.5)
hold on
plot(T,Y(:,1),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,1),'-r','LineWidth',1.5)
legend('Active (FT)','Active (IT)','Passive')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-SD','-djpeg','-r300')

fig = figure(4);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TFT,YFT(:,3),'-g','LineWidth',1.5)
hold on 
plot(T,Y(:,3),'-k','LineWidth',1.5)
plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
legend('Active (FT)','Active (IT)','Passive')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-TD','-djpeg','-r300')

fig = figure(5);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
subplot(2,2,1)
plot(TFT(1:end-1),K2(1:end-1,1),'-r','LineWidth',1.5)
hold on
plot(TFT,K(:,1)*ones(steps,1),'-g','LineWidth',1.5)
legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$K_1$','Interpreter','Latex','FontSize',12)
subplot(2,2,2)
plot(TFT(1:end-1),K2(1:end-1,2),'-r','LineWidth',1.5)
hold on
plot(TFT,K(:,2)*ones(steps,1),'-g','LineWidth',1.5)
legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$K_2$','Interpreter','Latex','FontSize',12)
subplot(2,2,3)
plot(TFT(1:end-1),K2(1:end-1,3),'-r','LineWidth',1.5)
hold on
plot(TFT,K(:,3)*ones(steps,1),'-g','LineWidth',1.5)
legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$K_3$','Interpreter','Latex','FontSize',12)
subplot(2,2,4)
plot(TFT(1:end-1),K2(1:end-1,4),'-r','LineWidth',1.5)
hold on
plot(TFT,K(:,4)*ones(steps,1),'-g','LineWidth',1.5)
legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$K_4$','Interpreter','Latex','FontSize',12)
% print('Active-TD-Dynamic-Costates','-djpeg','-r300')

%}

%% LQT


clear S;

tspan = [tf t0-stepsize];
S0 = zeros(16,1);
[TRiccati,S] = rk4fixed(@finiteLQRRiccati,tspan,S0,steps+1);

S = flipud(S);

%%
nu0T = [0; 0; 0; 0];
nuT = zeros(steps+1,4);
nuT(end,:) = nu0T';

tfiter = tf;

Qbar = Q - N*Rinv*N';
Abar = A - B*R*N';

% SVec = S(1,:);
% SMat = (reshape(SVec,[4,4]))';

%%
for i = steps + 1:-1:2
    i
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    t0iter = tfiter - stepsize;
    tspan = [tfiter t0iter];
    nu0Titer = nuT(i,:)';
    [~,nuTiter] = ode15s(@finiteLQTNu,tspan,nu0Titer);
    tfiter = t0iter;
    nuT(i-1,:) = nuTiter(end,:);
end

zs0T = 0;
zu0T = 0;
zsdot0T = 0;
zudot0T = 0;
zr0T = 0;

x10T = zs0T-zu0T;
x20T = zsdot0T;
x30T = zu0T-zr0T;
x40T = zudot0T;

x0T = [x10T; x20T; x30T; x40T];
YT = zeros(steps,4);
YT(1,:) = x0T';
TT = zeros(steps,1);
zacclT = zeros(steps,1);
K3 = zeros(steps,4);

t0iter = t0;

for i = 1:steps-1
    i
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    nuiter = nuT(i,:)';
    tfiter = t0iter + stepsize;
    tspan = [t0iter tfiter];
    x0Titer = YT(i,:)';
    [TTiter,YTiter] = ode15s(@car_lqt_finite,tspan,x0Titer);
    t0iter = tfiter;
    K3(i,:) = Rinv*(B'*SMat + N');
    [xdotT, zsddotT] = car_lqt_finite(TT(i),YT(i,:)');
    zacclT(i) = zsddotT;
    TT(i+1) = tfiter;
    YT(i+1,:) = YTiter(end,:);
end

% zuT = YT(:,3) + ZR;
% zsT = YT(:,1) + zuT;

%%
% fig = figure(6);
% set(fig,'Position',[1800 -320 1200 1000])
% clear title
% clear legend
% plot(TT,zsT,'-g','LineWidth',1.5)
% hold on
% plot(TPass(1:steps),zsPass(1:steps),'-r','LineWidth',1.5)
% plot(T,ZR,'-.b','LineWidth',1.5)
% title('Sprung Mass Deflection vs. Time')
% xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
% ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% legend('Active (Tracker)','Passive','Road Profile')
% set(legend,'Interpreter','Latex','FontSize',12)
% % print('Passive-SMD','-djpeg','-r300')
% 
% fig = figure(7);
% set(fig,'Position',[1800 -320 1200 1000])
% clear title
% clear legend
% plot(TT,zacclT,'-g','LineWidth',1.5)
% hold on
% plot(T,zaccl,'-k','LineWidth',1.5)
% plot(TPass(1:steps),zacclPass(1:steps),'-r','LineWidth',1.5)
% legend('Active (Tracker)','Passive')
% title('Sprung Mass Acceleration vs. Time')
% xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
% ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)
% % print('Passive-SMA','-djpeg','-r300')

fig = figure(8);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TT,YT(:,1),'-g','LineWidth',1.5)
hold on
% plot(T,Y(:,1),'-k','LineWidth',1.5)
% plot(TPass(1:steps),YPass(1:steps,1),'-r','LineWidth',1.5)
% legend('Active (Tracker)','Passive')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-SD','-djpeg','-r300')

fig = figure(9);
set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TT,YT(:,3),'-g','LineWidth',1.5)
hold on 
% plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
% legend('Active (Tracker)','Passive')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)

%}
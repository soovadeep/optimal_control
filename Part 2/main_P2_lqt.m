clear

global ks kt ms mu bs bt w Amp A B C L rho1 rho2 rho3 rho4 R N Q Rinv x10 x20 x30 x40 tf K SMat t0 nuiter Qbar Abar

kt = (704*10^3)/4; % N/m
ks = 15*10^3; % N/m
bs = 1400; % Ns/m  Check
bt = 0; % Ns/m
mu = 181/4; % kg
ms = 1814/4; % kg
rho1 = 1000000; % 0.4 
rho2 = 0.04; % 0.04
rho3 = 0.4; % 0.4
rho4 = 0.04; % 0.04
Amp = 0.05;
w = 0*2*pi;
t0 = 0;
tf = 1;
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

clear S;

tspan = [tf t0-stepsize];
S0 = zeros(16,1);
[TRiccati,S] = rk4fixed(@finiteLQRRiccati,tspan,S0,steps+1);

S = flipud(S);

%% LQT

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

% fig = figure(8);
% % set(fig,'Position',[1800 -320 1200 1000])
% clear title
% clear legend
% plot(TT,YT(:,1),'-g','LineWidth',1.5)
% hold on
% % plot(T,Y(:,1),'-k','LineWidth',1.5)
% % plot(TPass(1:steps),YPass(1:steps,1),'-r','LineWidth',1.5)
% % legend('Active (Tracker)','Passive')
% title('Suspension Deflection vs. Time')
% xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
% ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% % print('Passive-SD','-djpeg','-r300')
% 
% fig = figure(9);
% % set(fig,'Position',[1800 -320 1200 1000])
% clear title
% clear legend
% plot(TT,YT(:,3),'-g','LineWidth',1.5)
% hold on 
% % plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
% % legend('Active (Tracker)','Passive')
% title('Tire Deflection vs. Time')
% xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
% ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)

%}

%% Non-zero Setpoint

nu0NZ = [0; 0; 0; 0];
nuNZ = zeros(steps+1,4);
nuNZ(end,:) = nu0NZ';

tfiter = tf;

SVec = S(1,:);
SMat = (reshape(SVec,[4,4]))';

for i = steps + 1:-1:2
    i
%     SVec = S(i,:);
%     SMat = (reshape(SVec,[4,4]))';
    t0iter = tfiter - stepsize;
    tspan = [tfiter t0iter];
    nu0NZiter = nuNZ(i,:)';
    [~,nuNZiter] = ode15s(@finiteLQTNu,tspan,nu0NZiter);
    tfiter = t0iter;
    nuNZ(i-1,:) = nuNZiter(end,:);
end

zs0NZ = 0;
zu0NZ = 0;
zsdot0NZ = 0;
zudot0NZ = 0;
zr0NZ = 0;

x10NZ = zs0NZ-zu0NZ;
x20NZ = zsdot0NZ;
x30NZ = zu0NZ-zr0NZ;
x40NZ = zudot0NZ;

x0NZ = [x10NZ; x20NZ; x30NZ; x40NZ];
YNZ = zeros(steps,4);
YNZ(1,:) = x0NZ';
TNZ = zeros(steps,1);
zacclNZ = zeros(steps,1);
K4 = zeros(steps,4);

t0iter = t0;

for i = 1:steps-1
    i
%     SVec = S(i,:);
%     SMat = (reshape(SVec,[4,4]))';
    nuiter = nuNZ(i,:)';
    tfiter = t0iter + stepsize;
    tspan = [t0iter tfiter];
    x0NZiter = YNZ(i,:)';
    [~,YNZiter] = ode15s(@car_lqt_finite,tspan,x0NZiter);
    t0iter = tfiter;
    K4(i,:) = Rinv*(B'*SMat + N');
    [~, zsddotNZ] = car_lqt_finite(TNZ(i),YNZ(i,:)');
    zacclNZ(i) = zsddotNZ;
    TNZ(i+1) = tfiter;
    YNZ(i+1,:) = YNZiter(end,:);
end

% zuT = YT(:,3) + ZR;
% zsT = YT(:,1) + zuT;

%%
% fig = figure(14);
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
% fig = figure(15);
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

fig = figure(16);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TNZ,YNZ(:,1),'-g','LineWidth',1.5)
hold on
plot(TT,YT(:,1),'-k','LineWidth',1.5)
% plot(T,Y(:,1),'-k','LineWidth',1.5)
% plot(TPass(1:steps),YPass(1:steps,1),'-r','LineWidth',1.5)
legend('Active (NZ)','Active (Tracker)')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-SD','-djpeg','-r300')

fig = figure(17);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(TNZ,YNZ(:,3),'-g','LineWidth',1.5)
hold on 
plot(TT,YT(:,3),'-k','LineWidth',1.5)
% plot(TPass(1:steps),YPass(1:steps,3),'-r','LineWidth',1.5)
legend('Active (NZ)','Active (Tracker)')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)



%}
clear

global ks kt ms mu bs bt w Amp A B C L rho1 rho2 rho3 rho4 R N Q Rinv x10 x20 x30 x40 tf K SMat t0 nuiter Qbar Abar VVec PVec

load('Y_passive_static.mat')
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
tf = 10;
steps = tf*1000;
stepsize = (tf-t0)/steps;

A = [0 1 0 -1; -ks/ms -bs/ms 0 bs/ms; 0 0 0 1; ks/mu bs/mu -kt/mu -(bs+bt)/mu ];
B = [0;1/ms;0;-1/mu];
L = [0;0;-1;0];
% C = [1 0 0 0;0 1 0 0];
% C = eye(4);
C = [1 0 0 0];

R = 1/ms^2;
Rinv = 1/R;
N = [-ks/ms^2; -bs/ms^2; 0; bs/ms^2];
Q = [(ks^2/ms^2 + rho1)  bs*ks/ms^2            0      -bs*ks/ms^2;
     bs*ks/ms^2          (bs^2/ms^2 + rho2)    0      -bs^2/ms^2;
     0                   0                     rho3   0;
     -bs*ks/ms^2         -bs^2/ms^2            0      (bs^2/ms^2 + rho2)];
Qbar = Q - N*Rinv*N';
Abar = A - B*R*N';

%% Final States Fixed 

tspan = [tf t0-stepsize];
S0 = zeros(16,1);
[~, S] = rk4fixed(@finalStateFixed_S, tspan, S0, steps+1);

S = flipud(S);

%%
% V0T = [0; 0; 0; 0];
V0 = C';
V = zeros(steps+1,4);
V(end,:) = V0';

tfiter = tf;

for i = steps + 1:-1:2
    i
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    t0iter = tfiter - stepsize;
    tspan = [tfiter t0iter];
    V0iter = V(i,:)';
    [~,Viter] = ode23s(@finalStateFixed_V,tspan,V0iter);
    tfiter = t0iter;
    V(i-1,:) = Viter(end,:);
end

%%

P0 = [0];
% P0 = C';
P = zeros(steps+1,1);
P(end,:) = P0';

tfiter = tf;

for i = steps + 1:-1:2
    i
%     SVec = S(i,:);
%     SMat = (reshape(SVec,[4,4]))';
    VVec = V(i,:)';
    t0iter = tfiter - stepsize;
    tspan = [tfiter t0iter];
    P0iter = P(i)';
    [~,Piter] = ode15s(@finalStateFixed_P,tspan,P0iter);
    tfiter = t0iter;
    P(i-1) = Piter(end);
end

%%
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

Y = zeros(steps,4);
T = zeros(steps,1);
zaccl = zeros(steps,1);
force = zeros(steps,1);
Y(1,:) = x0';
K2 = zeros(steps,4);

t0iter = t0;

for i = 1:steps-1
    i
    SVec = S(i,:);
    SMat = (reshape(SVec,[4,4]))';
    VVec = V(i,:)';
    PVec = P(i);
    tfiter = t0iter + stepsize;
    tspan = [t0iter tfiter];
    x0iter = Y(i,:)';
    [Titer,Yiter] = ode15s(@finalStateFixed,tspan,x0iter); %400
    t0iter = tfiter;
    K2(i,:) = Rinv*(B'*SMat + N');
    [Ydot, zsddot, u] = finalStateFixed(T(i),Y(i,:)');
    zaccl(i) = zsddot;
    force(i) = u;
    T(i+1) = tfiter;
    Y(i+1,:) = Yiter(end,:);
end

ZR = Amp*sin(w*T);

zu = Y(:,3) + ZR;
zs = Y(:,1) + zu;

steps = tf*1000;

%%
fig = figure(1);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,zs,'-g','LineWidth',1.5)
hold on
plot(TPass(1:steps),zsPass(1:steps),'-.r','LineWidth',1.5)
plot(T,ZR,'-.b','LineWidth',1.5)
title('Sprung Mass Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
legend('Active (FSF)','Passive','Road Profile')
% set(legend,'Interpreter','Latex','FontSize',12)
% print('Passive-SMD-FSF','-djpeg','-r300')

fig = figure(2);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,zaccl,'-g','LineWidth',1.5)
hold on
plot(TPass(1:steps),zacclPass(1:steps),'-.r','LineWidth',1.5)
legend('Active (FSF)','Passive')
title('Sprung Mass Acceleration vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$\ddot{Z}_s\hspace{0.05in}(m/s^2)$','Interpreter','Latex','FontSize',12)
% print('Passive-SMA-FSF','-djpeg','-r300')

fig = figure(3);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,Y(:,1),'-g','LineWidth',1.5)
hold on
plot(TPass(1:steps),YPass(1:steps,1),'-.r','LineWidth',1.5)
legend('Active (FSF)','Passive')
title('Suspension Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_s - Z_u\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)
% print('Passive-SD-FSF','-djpeg','-r300')

fig = figure(4);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,Y(:,3),'-g','LineWidth',1.5)
hold on 
plot(TPass(1:steps),YPass(1:steps,3),'-.r','LineWidth',1.5)
legend('Active (FSF)','Passive')
title('Tire Deflection vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$Z_u - Z_r\hspace{0.05in}(m)$','Interpreter','Latex','FontSize',12)

%%
fig = figure(5);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
subplot(2,2,1)
plot(T,V(1:end-1,1),'-r','LineWidth',1.5)
hold on
% plot(TT,nuiter(1)*ones(steps,1),'-g','LineWidth',1.5)
% legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$V_1$','Interpreter','Latex','FontSize',12)
subplot(2,2,2)
plot(T,V(1:end-1,2),'-r','LineWidth',1.5)
hold on
% plot(TT,nuiter(2)*ones(steps,1),'-g','LineWidth',1.5)
% legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$V_2$','Interpreter','Latex','FontSize',12)
subplot(2,2,3)
plot(T,V(1:end-1,3),'-r','LineWidth',1.5)
hold on
% plot(TT,nuiter(3)*ones(steps,1),'-g','LineWidth',1.5)
% legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$V_3$','Interpreter','Latex','FontSize',12)
subplot(2,2,4)
plot(T,V(1:end-1,4),'-r','LineWidth',1.5)
hold on
% plot(TT,nuiter(4)*ones(steps,1),'-g','LineWidth',1.5)
% legend('Finite Time','Infinite Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
ylabel('$V_4$','Interpreter','Latex','FontSize',12)
% print('V-FSF','-djpeg','-r300')
%}

fig = figure(6);
% set(fig,'Position',[1800 -320 1200 1000])
clear title
clear legend
plot(T,P(1:end-1),'-g','LineWidth',1.5)
hold on 
% plot(TPass(1:steps),YPass(1:steps,3),'-.r','LineWidth',1.5)
% legend('Active (FSF)','Passive')
title('P vs. Time')
xlabel('$Time\hspace{0.05in}(s)$','Interpreter','Latex','FontSize',12)
% ylabel('$P$','Interpreter','Latex','FontSize',12)
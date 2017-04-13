%Simulating function
clear
clc
global ks kt ms mu bs bt zr w Amp

Amp=0.5;
w=0.5;

t0=0;
tf=100;

zs0=0;
zu0=0;
zsdot0=0;
zudot0=0;
zr0=0;

x10=zs0-zu0;
x20=zsdot0;
x30=zu0-zr0;
x40=zudot0;

x0=[x10;x20;x30;x40];

tspan=[t0 tf];
[T,Y]=rk4fixed(@car,tspan,x0,10000);

u=0;
ZR=Amp*sin(w*T);

zu=Y(:,3)+ZR;
zs=Y(:,1)+zu;
figure(1)
plot(T,zs,'r',T,ZR,'b');
title('Zs vs T');
figure(2)
plot(T,zu);
title('Zu vs T');

figure(3)
plot(T,ZR);
title('ZR vs T');

figure(4)
plot(T,Y(:,1));

%x1=zs-zu
%x2=zsdot
%x3=zu-zr
%x4=zudot
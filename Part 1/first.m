clear
clc

syms ks ms bs ms bs bt kt mu
A = [0 1 0 -1; -ks/ms -bs/ms 0 bs/ms; 0 0 0 1; ks/mu bs/mu -kt/mu -(bs)/mu ];
% B = [0;1/ms;0;-1/mu];
% L = [0;0;-1;0];
% C = [1 0 -1 0;0 1 0 0];
charpoly(A)
clear,clc,close all;
xmin=0;
xmax=1;
tmin=0;
tmax=10;
Nx=256;
Nt=4096;
X=linspace(xmin,xmax,Nx);
T=linspace(tmin,tmax,Nt);
[Xm,Tm]=meshgrid(X,T);
M=150;
x0_vals = linspace(xmin,xmax,M);
t0_vals = linspace(tmin,tmax,M);
tau = 0.5;
tcenter = 3;
A=2;

dx = X(2)-X(1);
dt = T(2)-T(1);

save("datas.mat")
disp("Initialized!")
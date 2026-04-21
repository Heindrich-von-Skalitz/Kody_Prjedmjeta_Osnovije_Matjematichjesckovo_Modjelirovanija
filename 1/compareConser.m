clear,clc,close all;
cons=load("conser.mat","u");
noncons=load("nonconser.mat","u");
load("datas.mat")
u1=cons.u;
u2=noncons.u;
Delta = zeros(Nx,Nt);
for i=1:Nx
    for j=1:Nt
        Delta(i,j)=abs(u1(i,j)-u2(i,j));
    end
end

figure;
hold on; grid on; box on;
contourf(X, T, Delta', 20, 'LineColor', 'none');
xlabel('x');
ylabel('t');
title('Δu(x,t)');
colorbar;
colormap(jet);
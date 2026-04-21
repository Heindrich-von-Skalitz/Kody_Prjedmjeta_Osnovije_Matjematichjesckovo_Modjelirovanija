clear,clc,close all;
load("datas.mat")
u=zeros(Nx,Nt);
u(:,1)=X;

%Гауссовый пакет
for j = 1:Nt
    u(1,j)=A*exp(-(T(j)-tcenter)^2/tau^2);
end

%% Кор Код
for j=1:(Nt-1)
    for i = 2:Nx
        u(i,j+1) = u(i,j) - dt/(1+u(i,j)^2)*(u(i,j) - u(i-1,j))/dx;
    end
end

%% Графики
figure;
hold on; grid on; box on;
contourf(X, T, u', 20, 'LineColor', 'none');
xlabel('x');
ylabel('t');
title('u(x,t) с характеристики');
plot(NaN, NaN, 'b-', 'LineWidth', 1, 'DisplayName', '(t*=0)');
plot(NaN, NaN, 'r-', 'LineWidth', 1, 'DisplayName', '(x*=0)');
xlim([xmin,xmax]);
ylim([tmin,tmax]);
colorbar;
colormap(jet);

for i = 1:length(x0_vals)
    x0 = x0_vals(i);
    t = linspace(tmin, tmax, Nt);
    x = x0 + t./(1+x0^2);
    plot(x, t, 'b-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
for j = 1:length(t0_vals)
    t0 = t0_vals(j);
    t = linspace(t0, tmax, Nt);
    x=(t-t0)/(1+A^2*exp(-2*(t0-tcenter)^2/tau^2));
    plot(x, t, 'w-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

figure;%Без характеристики
hold on; grid on; box on;
contourf(X, T, u', 20, 'LineColor', 'none');
xlabel('x');
ylabel('t');
title('u(x,t)');
xlim([xmin,xmax]);
ylim([tmin,tmax]);
colorbar;
colormap(jet);

figure;
plot3(Xm,Tm,u)

save("nonconser.mat","u")
clear,clc,close all;
load("datas.mat")
u=zeros(Nx,Nt);
u(:,1)=X;

epsilon = 1E-10;
maxTimes = 1E4;

%Гауссовый пакет
for j = 1:Nt
    u(1,j)=A*exp(-(T(j)-tcenter)^2/tau^2);
end

inp = 0;
for n=1:Nt-1
    inp = inp+u(1,n)*dt;
end

%% Кор Код
for j=1:(Nt-1)
    for i = 2:Nx
        u(i,j+1) = u(i,j) - dt*(atan(u(i,j)) - atan(u(i-1,j)))/dx;
    end
end
% for j=1:(Nt-1)
%     for i = 2:Nx
%         v = u(i,j);
%         for k = 1:maxTimes
%             F = v - u(i,j) + dt*(atan(v) - atan(u(i-1,j+1)))/dx;
%             dF = 1+dt/(dx*(1+v^2));
%             v = v-F/dF;
%             if abs(u(i,j+1)-v)<epsilon
%                 u(i,j+1) = v;
%                 break
%             end
%             u(i,j+1) = v;
%         end
%     end
% end

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

save("gauss.mat","u")

balan = 0;
for m = 1:Nx
    balan = balan + dx * (u(m,Nt) - u(m,1));
end
for n = 1:Nt-1
    balan = balan + dt * atan(u(Nx,n));
end

[ux, ut] = gradient(u, dx, dt);
figure;
hold on; grid on; box on;
contourf(X, T, ut', 20, 'LineColor', 'none');
xlabel('x');
ylabel('t');
title('u_t(x,t)');
xlim([xmin,xmax]);
ylim([tmin,tmax]);
colorbar;
colormap(jet);

figure;
hold on; grid on; box on;
contourf(X, T, ux', 20, 'LineColor', 'none');
xlabel('x');
ylabel('t');
title('u_x(x,t)');
xlim([xmin,xmax]);
ylim([tmin,tmax]);
colorbar;
colormap(jet);
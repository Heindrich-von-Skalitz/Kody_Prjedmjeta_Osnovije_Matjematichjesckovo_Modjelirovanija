clear,clc,close all;
load("datas.mat")
load("newton.mat")
balan=0;
dudt=zeros(Nx,Nt);
dudx=zeros(Nx,Nt);
balanc=zeros(Nx,Nt);

% for m=1:Nx
%     for n=1:Nt
%         if and(not(n==1),not(n==Nt))
%             dudt(m,n)=0.5*(u(m,n+1)-u(m,n-1))/dt;
%         elseif(n==1)
%             dudt(m,n)=(u(m,n+1)-u(m,n))/dt;
%         elseif(n==Nt)
%             dudt(m,n)=(u(m,n)-u(m,n-1))/dt;
%         end
%         if and(not(m==1),not(m==Nx))
%             dudt(m,n)=0.5/(1+u(m,n))*(u(m+1,n)-u(m-1,n))/dt;
%         elseif(m==1)
%             dudt(m,n)=1/(1+u(m,n))*(u(m+1,n)-u(m,n))/dt;
%         elseif(m==Nx)
%             dudt(m,n)=1/(1+u(m,n))*(u(m,n)-u(m-1,n))/dt;
%         end
%     end
% end
% for m=1:Nx
%     for n=1:Nt
%         balanc(m,n)=dx*dt*(dudt(m,n)+dudx(m,n));
%         balan=balan+balanc(m,n);
%     end
% end
% figure;
% hold on; grid on; box on;
% contourf(X, T, balanc', 20, 'LineColor', 'none');
% xlabel('x');
% ylabel('t');
% title('Result of balance');
% colorbar;
% colormap(jet);
balan = 0;
for m = 1:Nx
    balan = balan + dx * (u(m,Nt) - u(m,1));
end
for n = 1:Nt-1
    balan = balan + dt * atan(u(Nx,n));
end
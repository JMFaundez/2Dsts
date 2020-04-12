close all
clear all
addpath matlab_script/

% Load an unroll data
L = load('stsINT6.mat');
data = L.data_int;

x = data(:,:,1);
y = data(:,:,2);
U = data(:,:,3);
V = data(:,:,4);
W = data(:,:,5);
uu = data(:,:,7);
vv = data(:,:,8);
ww = data(:,:,9);
uv = data(:,:,11);
[ny,nx] = size(x);

file_der = 'dertan_der0.f00001';
nelx = 300;
nely = 22;
[der_data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(file_der);
[x_d,y_d,dtdn,dwdn] = reshapenek(der_data,nelx,nely);

alpha = zeros(ny,nx);
dx = x(2,:) - x(1,:);
dy = y(2,:) - y(1,:);
L = sqrt(dx.^2 + dy.^2);
a = acos(dy./L);
s = zeros(nx,1);
s(2:end,1) = cumsum(sqrt((x(1,2:end)-x(1,1:end-1)).^2+(y(1,2:end)-y(1,1:end-1)).^2));
for i=1:nx
  alpha(:,i) = a(i);
end

cos_a = cos(alpha);
sin_a = sin(alpha);

% Create normal and tangential coordinates
S = zeros(ny,nx);
N = zeros(ny,nx);
for n=1:ny
  S(n,:) = x(1,:);
  dx = x(n,:) - x(1,:);
  dy = y(n,:) - y(1,:);
  dh = sqrt(dx.^2+dy.^2);
  N(n,:) = dh;
end

% Compute tangential components
tt = uu.*cos_a.^2 + vv.*sin_a.^2 + 2*uv.*cos_a.*sin_a;
T = U.*cos_a + V.*sin_a;
pp = tt -T.^2;
pp(abs(pp)<1e-14)=0;
p = sqrt(pp);

wwp = ww-W.^2;
wwp(abs(wwp)<1e-14)=0;
wp = sqrt(wwp);

% Compute cf
p_max = zeros(nx,1);
n_max = zeros(nx,1);
for i=1:nx
  [p_max(i), jmax] = max(p(:,i));
  n_max(i) = N(jmax,i);
end

in0 = 490;
x(1,in0)
Re = 5.33333e5;
s = s - s(in0);
s = s*Re;
cf1 = zeros(nx,1);
cf2 = zeros(nx,1);
fin = 10;
for i=1:nx
  cf1(i) = 2*(-p(1,i) + p(fin,i))/(N(fin,i) - N(1,i))/Re;
  cf2(i) = 2*(-wp(1,i) + wp(fin,i))/(N(fin,i) - N(1,i))/Re;
end
%cf1 = smoothdata(cf1,'Gaussian',10);

figure()
subplot(222)
contourf(S(:,in0:end),N(:,in0:end)*sqrt(Re),p(:,in0:end), 'LineStyle','none')
xlim([0,0.4])
ylim([0,0.01*sqrt(Re)])
ylabel('$n/\delta$','Interpreter','latex')
xlabel('$x/c$','Interpreter','latex')
title("$\sqrt{\overline{u_{t}'u_{t}'}}$",'Interpreter','latex')
colorbar('manual','Position',[0.93, 0.73, 0.02, 0.20])

subplot(221)
plot(S(1,in0:end),p_max(in0:end),'linewidth',1.5)
xlim([0,0.4])
grid on
ylabel("$\left(\sqrt{\overline{u_{t}'u_{t}'}}\right)_{max}$",...
       'Interpreter','latex', 'FontSize',20)

stp = 10;
subplot(223)
hold on
plot(S(1,in0:stp:end),cf1(in0:stp:end),'linewidth',1.5)
plot(S(1,in0:stp:end),2/Re*dtdn(1,in0:stp:end),'linewidth',1.5)
xlabel("$x/c$",'Interpreter','latex','FontSize',14)
xlim([0,0.4])
grid on
ylabel("$c_{fx}$",'Interpreter','latex', 'FontSize',20)

subplot(224)
plot(S(1,in0:end),cf2(in0:end),'linewidth',1.5)
grid on
ylabel('$c_{fz}$','Interpreter','latex', 'FontSize',20)
xlabel("$x/c$",'Interpreter','latex','FontSize',16)
xlim([0,0.4])


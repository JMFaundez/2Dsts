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


alpha = zeros(ny,nx);
dx = x(2,:) - x(1,:);
dy = y(2,:) - y(1,:);
L = sqrt(dx.^2 + dy.^2);
a = acos(dy./L);
for i=1:nx
  alpha(:,i) = a(i);
end

cos_a = cos(alpha);
sin_a = sin(alpha);

% Compute tangential components
tt = uu.*cos_a.^2 + vv.*sin_a.^2 + 2*uv.*cos_a.*sin_a;
T = U.*cos_a + V.*sin_a;
pp = tt -T.^2;
pp(abs(pp)<1e-14)=0;
p = sqrt(pp);

wwp = ww-W.^2;
wwp(abs(wwp)<1e-14)=0;
wp = sqrt(wwp);


to_mesh = zeros(ny,nx,6);
to_mesh(:,:,1) = x;
to_mesh(:,:,2) = y;
to_mesh(:,:,3) = p;
to_mesh(:,:,4) = wp;

to_nek = demeshnek(to_mesh,10);

f_name = "pert0.f00001";
[dum_data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('fringe_m90.f00008');


emode = 'le';
status = writenek(f_name,to_nek,lr1,elmap,time,istep,fields,emode,wdsz,etag) 


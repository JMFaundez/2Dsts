clear all

set = 4;
filename = ['stsINT',num2str(set),'.mat'];
L = load(filename);
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

data_int = zeros(ny,nx,5);

upert = sqrt((abs(uu-U.^2) + abs(vv-V.^2) + abs(ww-W.^2))/3);
U_av = sqrt(U.^2 + V.^2 + W.^2);
TI_1 = upert;
TI_2 = upert./U_av;

data_int(:,:,1) = x;
data_int(:,:,2) = y;
data_int(:,:,4) = TI_1;
data_int(:,:,5) = TI_2;
fileoutput = ['intTI',num2str(set),'_5p.mat'];
save(fileoutput, 'data_int')

figure()
subplot(211)
contourf(x,y,TI_1,'LineColor','none')
colorbar()

in = 1;
subplot(212)
contourf(x(in:end,:),y(in:end,:),TI_2(in:end,:),'LineColor','none')
caxis([0 0.12])
colorbar()


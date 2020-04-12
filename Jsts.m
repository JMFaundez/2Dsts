close all
clear all
clc
addpath matlab_script/

%% Load fils
% this is the file with the right gll ordering
input_mesh='./fringe_m90.f00008';

% Number of elements x and y direction
nelx = 300;
nely = 22;

% Choose the set of files to average and interpolate
set = 6;
[files, T] = set_files(set);


[data_mesh,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input_mesh);
[xx,yy] = reshapenek(data_mesh,nelx,nely);
[nx,ny] = size(xx);

%% Interpolate and average data
N = length(files);
data_int = zeros(nx,ny,46);
tT = 0;
tic
for i=1:N
  i
  input_sts = files(i);
  ti = T(i,:);
  [data_sts,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input_sts);
  data_i = interpolate_data(xx,yy,data_sts,ti);
  tT = tT +  ti(2) - ti(1);
  data_int = data_int + data_i;
end
data_int = data_int/tT;
data_int(:,:,1) = xx;
data_int(:,:,2) = yy;
toc

output_sts = ['stsINT',num2str(set),'.mat']

save(output_sts,'data_int')

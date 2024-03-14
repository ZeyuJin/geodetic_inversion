%experimentally read the configuration file
%11/3/2022, Xiaoyu Zou, x3zou@ucsd.edu
clear
clc
% addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/sign_mask')
% addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/ASC100/LOS2')
% addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/DES5/LOS2')
% setenv('PATH',[getenv('PATH'),':/opt/homebrew/bin']);  % add the path of GMT
% configtest='configtest.txt';
% fid = fopen(configtest);
% tmp_txt = fgetl(fid);
% files=zeros(1,1);
% files=string(files);
% while tmp_txt ~=-1
%     tmp_txt = fgetl(fid);
%     skip=find(tmp_txt=='#');
%     if ~isempty(skip)
%         continue
%     end
%     files=[files tmp_txt];
% end
% files(1)=[];
% files(end)=[];
% fclose(fid)


clear
clc
addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/sign_mask')
addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/ASC100/LOS2')
addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/DES5/LOS2')
setenv('PATH',[getenv('PATH'),':/opt/homebrew/bin']);  % add the path of GMT
configtest='configtest2.txt';
fid = fopen(configtest);
tmp_txt = fgetl(fid);
para=zeros(1,1);
para=string(para);
while tmp_txt ~=-1
    tmp_txt = fgetl(fid);
    skip=find(tmp_txt=='#');
    if ~isempty(skip)
        continue
    end
    para=[para tmp_txt];
end
para(1)=[];
para(end)=[];
para=str2double(para);
fclose(fid)
lonf=para(1);
latf=para(2);
ref_lon=para(3);
threshold=para(4);
lonc=para(5);
latc=para(6);
Nmin=para(7);
Nmax=para(8);
dip_change_id=para(9):para(10);
dip_angle=[para(11) para(12) para(13) para(14) para(15)];
iter_step=para(16);
iter_step2=para(17);
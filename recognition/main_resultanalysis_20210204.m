clc
clear all
close all

datapath='./output/';
L=dir([datapath,'/*.txt']);
for ii=1:length(L)
    filename=L(ii).name;
    filepath=[datapath,filename];
    
    fid=fopen(filepath,'r');
    tline=fgetl(fid);
    ind=strfind(tline,'accuracy=');
    accu(ii)=str2double(tline(ind+length('accuracy='):end));

    fclose(fid);
end

mean(accu)

find(accu<0.9)
accu(find(accu<0.9))

mean(accu(accu>0.9))

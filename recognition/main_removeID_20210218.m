% remove manual annotated cell name - for test
% by Lei Qu 20210218

clc
clear all
close all

datapath='./input_withID/';

L=dir([datapath,'/*.ano.ano.txt']);
for iter_file=1:length(L)
    filename_apo=L(iter_file).name;
    filepath_apo=[datapath,filename_apo];

    apo=load_v3d_pointcloud_file(filepath_apo);
    fprintf('load apo file:[%s]\n',filename_apo);
    fprintf('\t[%d] cells load\n',length(apo));

    %remove manual ID 
    for i=1:length(apo)
        cellname=apo{i}.name;
        if ~(strcmpi(cellname,'f1') || strcmpi(cellname,'f2') || strcmpi(cellname,'f3') || strcmpi(cellname,'f4') || strcmpi(cellname,'f5') || ...
                strcmpi(cellname,'f6') || strcmpi(cellname,'f7') || strcmpi(cellname,'f8') || strcmpi(cellname,'f9') || strcmpi(cellname,'f10') || ...
                contains(cellname,'nouse') || contains(cellname,'NOUSE'))
            apo{i}.name='';
        end
    end

    %save 
    output_filename=['./input_noID/',filename_apo];
    fp = fopen(output_filename, 'w');
    for i=1:length(apo)
        S=apo{i};
        fprintf(fp, '%d,%s,%s,%s,%d,%d,%d,%.2f,%.2f,%.2f,%d,%d\n', ...
        S.n, ...
        strtrim(S.orderinfo), ...
        strtrim(S.name), ...
        strtrim(S.comment), ...
        S.z, ...
        S.x, ...
        S.y, ...
        S.pixmax, ...
        S.intensity, ...
        S.sdev, ...
        S.volsize, ...
        S.mass);           
    end
    fclose(fp);
    fprintf('Save predict result to [%s] done.\n', output_filename);
 
    
end


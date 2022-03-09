% C.elegant cell recognition based on manual segmentation
% by Lei Qu 20200315
% modified by Lei Qu @ 20210202
%       1. predict cell's ID (cell num = 558)
%       2. given manual recog result, calculate recognition accuracy
% modified by Lei Qu @ 20210217
%       1. deal cell num < 558
%       2. output cell recog fidelity -> add new term apo.fidelity
%       3. preserve some cell's manual recog result(comment=*FIX*)

function main()
    clc
    clear all
    close all
    addpath('./RPM/');
    addpath('./matlab_io_basicdatatype/');
    addpath('./bipartite');
    

    filepath_apo_sub_gt="D:\Vaa3d\vaa3d_tools\released_plugins\v3d_plugins\C_elegent\c_elegent2\main_c_elegent\itera1\C07E3.6_SD1794_20140116_0001\affine_marker_label.txt";
    apo_sub_gt=load_v3d_apo_file(filepath_apo_sub_gt);
    filepath_apo_atlas='.\traindata\b_atlas_train100pwaff_20200417.apo';
    apo_atlas=load_v3d_pointcloud_file(filepath_apo_atlas);
    X_tar=-1*ones(10,3);
    X_sub=-1*ones(10,3);

    for i=1:10
        X_tar(i,1)=apo_atlas{i}.x;
        X_tar(i,2)=apo_atlas{i}.y;
        X_tar(i,3)=apo_atlas{i}.z;
        X_sub(i,1)=apo_sub_gt{i}.x;
        X_sub(i,2)=apo_sub_gt{i}.y;
        X_sub(i,3)=apo_sub_gt{i}.z;
    end
    lamda1=0.1;
      [w,d,k]  = update_tps( X_sub, X_tar,lamda1);

      vx = update_marker(X_sub,X_sub,w,d);
      
    fid = fopen(".\output\result_marker.txt", 'w');
    for i=1:length(10)
    fprintf(fid, ' %5.3f, %5.3f, %5.3f,0,0,,,255,255,255\n', vx(i,1), vx(i,2), vx(i,3));    
    end
    fclose(fid);
    

end
        



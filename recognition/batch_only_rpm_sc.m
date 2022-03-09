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
addpath('./RPM/codegen/mex/');

%load trained cell pos std
load('.\traindata\man2atlas_pwaffine_meanstd_norm.mat','arr_pos_std_pwaff');
arr_pos_var_pwaff=arr_pos_std_pwaff.^2; clear arr_pos_std_pwaff

load('.\traindata\man2atlas_pwaffine_meanstd_norm_sc.mat','arr_pos_mean_pwaff');%arr_pos_mean_pwaff,arr_pos_std_pwaff
arr_pos_var_pwaff_sc=arr_pos_mean_pwaff.^2; clear arr_pos_mean_pwaff

%'noID'(output recog result no accu), 'withID'(output recog accu)
datapath='Y:\2.postgrad\16-liyuanyuan\C_elegent\data\wangmeng_pre_mask\segmentation2point_size_shape_point_all\0\';
data_gt="Y:\2.postgrad\16-liyuanyuan\C_elegent\data\cellseganno_4train_20200330\apo_annotation_index_new\";
datapath_file="Y:\2.postgrad\16-liyuanyuan\C_elegent\data\wangmeng_pre_mask\segmentation2point_size_shape_point_all\";
L=dir([datapath,'/*.apo']);
for iter_file=1:length(L)
    if(iter_file==81)
        continue;
    end
    filename_sub=L(iter_file).name;
    filepath_apo_sub=[datapath,filename_sub];
    filepath_apo_sub_gt=strcat(data_gt,filename_sub(1:length(filename_sub)-8));
    filepath_apo_sub_gt=strcat(filepath_apo_sub_gt,'_mask.apo');
    apo_sub_gt=load_v3d_apo_file(filepath_apo_sub_gt);
    
    
    fprintf('==========================================================\n');
    fprintf('#(%d):%s\n',iter_file,filename_sub);
    fprintf('==========================================================\n');
    
    fprintf('[1] Load data ...\n');
    % Load sub apo file
    apo_sub=load_v3d_apo_file(filepath_apo_sub);
    file_name=strcat(datapath_file,"\1\",filename_sub);
    apo_sub1=load_v3d_apo_file(file_name);
    file_name=strcat(datapath_file,"\2\",filename_sub);
    apo_sub2=load_v3d_apo_file(file_name);
    file_name=strcat(datapath_file,"\3\",filename_sub);
    apo_sub3=load_v3d_apo_file(file_name);
    file_name=strcat(datapath_file,"\4\",filename_sub);
    apo_sub4=load_v3d_apo_file(file_name);
    file_name=strcat(datapath_file,"\5\",filename_sub);
    apo_sub5=load_v3d_apo_file(file_name);
    file_name=strcat(datapath_file,"\6\",filename_sub);
    apo_sub6=load_v3d_apo_file(file_name);
    fprintf('\tload sub apo file:[%s]\n',filename_sub);
    fprintf('\t\t[%d] cells load\n',length(apo_sub));
    clear ind filepath_apo_sub
    %Load atlas apo file
    filepath_apo_atlas='.\traindata\atlas_with_size_intensit.apo';
    apo_atlas=load_v3d_pointcloud_file(filepath_apo_atlas);
    apo_atlas1=load_v3d_pointcloud_file('.\traindata\update_atlas_size1.apo');
    apo_atlas2=load_v3d_pointcloud_file('.\traindata\update_atlas_size2.apo');
    apo_atlas3=load_v3d_pointcloud_file('.\traindata\update_atlas_size3.apo');
    apo_atlas4=load_v3d_pointcloud_file('.\traindata\update_atlas_size4.apo');
    apo_atlas5=load_v3d_pointcloud_file('.\traindata\update_atlas_size5.apo');
    apo_atlas6=load_v3d_pointcloud_file('.\traindata\update_atlas_size6.apo');
    ind=strfind(filepath_apo_atlas,'\');
    filename_atlas=filepath_apo_atlas(ind(end)+1:end);
    fprintf('\tload atlas apo file:[%s]\n',filename_atlas);
    fprintf('\t\t[%d] cells load\n',length(apo_atlas));
    clear ind filepath_apo_atlas filename_atlas
    
    fprintf('[4] Do cell recog...\n');
    nvalidcell=length(apo_sub);
    % reformat valid apo data to arr
    X_sub_gt=-1*ones(3,558);
    X_tar=-1*ones(3,558);
    X_tar_size=-1*ones(3,558);
    X_tar_intensity=-1*ones(1,558);
    
    %     X_sub=-1*ones(3,length(apo_sub));
    %     X_sub_size=-1*ones(3,length(apo_sub));
    %     X_sub_intensity=-1*ones(1,length(apo_sub));
    %     X_sub_shape=-1*ones(3,length(apo_sub));
    for i=1:6
        X_tar_shape_all{i}=-1*ones(3,558);
        X_sub_shape_all{i}=-1*ones(3,length(apo_sub));
    end
    
    for i=1:558
        X_tar(1,i)=apo_atlas{i}.x;
        X_tar(2,i)=apo_atlas{i}.y;
        X_tar(3,i)=apo_atlas{i}.z;
        X_tar_size(1,i)=apo_atlas{i}.pixmax;
        X_tar_size(2,i)=apo_atlas{i}.mass;
        X_tar_size(3,i)=apo_atlas{i}.sdev;
        X_tar_intensity(1,i)=apo_atlas{i}.volsize;
        
        X_sub_gt(1,i)=apo_sub_gt{i}.x;
        X_sub_gt(2,i)=apo_sub_gt{i}.y;
        X_sub_gt(3,i)=apo_sub_gt{i}.z;
        
         X_tar_shape_all{1}(1,i)=apo_atlas1{i}.x;
         X_tar_shape_all{1}(2,i)=apo_atlas1{i}.y;
         X_tar_shape_all{1}(3,i)=apo_atlas1{i}.z;
       
        X_tar_shape_all{2}(1,i)=apo_atlas2{i}.x;
        X_tar_shape_all{2}(2,i)=apo_atlas2{i}.y;
        X_tar_shape_all{2}(3,i)=apo_atlas2{i}.z;
     
         X_tar_shape_all{3}(1,i)=apo_atlas3{i}.x;
         X_tar_shape_all{3}(2,i)=apo_atlas3{i}.y;
         X_tar_shape_all{3}(3,i)=apo_atlas3{i}.z;

        X_tar_shape_all{4}(1,i)=apo_atlas4{i}.x;
        X_tar_shape_all{4}(2,i)=apo_atlas4{i}.y;
        X_tar_shape_all{4}(3,i)=apo_atlas4{i}.z;

         X_tar_shape_all{5}(1,i)=apo_atlas5{i}.x;
         X_tar_shape_all{5}(2,i)=apo_atlas5{i}.y;
         X_tar_shape_all{5}(3,i)=apo_atlas5{i}.z;

        X_tar_shape_all{6}(1,i)=apo_atlas6{i}.x;
        X_tar_shape_all{6}(2,i)=apo_atlas6{i}.y;
        X_tar_shape_all{6}(3,i)=apo_atlas6{i}.z;

    end
    n=0;
    for i=1:length(apo_sub)
        %         if((abs(apo_sub1{i}.x-apo_sub2{i}.x)<5&&abs(apo_sub3{i}.y-apo_sub4{i}.y)<5&&abs(apo_sub5{i}.z-apo_sub6{i}.z)<5))
        if(0)
            continue;
        else
            n=n+1;
            index_seg(n)=apo_sub{i}.n;
            X_sub(1,n)=apo_sub{i}.x;
            X_sub(2,n)=apo_sub{i}.y;
            X_sub(3,n)=apo_sub{i}.z;
            X_sub_size(1,n)=apo_sub{i}.pixmax;
            X_sub_size(2,n)=apo_sub{i}.mass;
            X_sub_size(3,n)=apo_sub{i}.sdev;
            X_sub_intensity(1,n)=apo_sub{i}.volsize;
            
            X_sub_shape_all{1}(1,n)=apo_sub1{i}.x;
            X_sub_shape_all{1}(2,n)=apo_sub1{i}.y;
            X_sub_shape_all{1}(3,n)=apo_sub1{i}.z;

             X_sub_shape_all{2}(1,n)=apo_sub2{i}.x;
             X_sub_shape_all{2}(2,n)=apo_sub2{i}.y;
             X_sub_shape_all{2}(3,n)=apo_sub2{i}.z;

            X_sub_shape_all{3}(1,n)=apo_sub3{i}.x;
            X_sub_shape_all{3}(2,n)=apo_sub3{i}.y;
            X_sub_shape_all{3}(3,n)=apo_sub3{i}.z;

            X_sub_shape_all{4}(1,n)=apo_sub4{i}.x;
            X_sub_shape_all{4}(2,n)=apo_sub4{i}.y;
            X_sub_shape_all{4}(3,n)=apo_sub4{i}.z;

            X_sub_shape_all{5}(1,n)=apo_sub5{i}.x;
            X_sub_shape_all{5}(2,n)=apo_sub5{i}.y;
            X_sub_shape_all{5}(3,n)=apo_sub5{i}.z;

            X_sub_shape_all{6}(1,n)=apo_sub6{i}.x;
            X_sub_shape_all{6}(2,n)=apo_sub6{i}.y;
            X_sub_shape_all{6}(3,n)=apo_sub6{i}.z;
        end
        
    end
    
    
    % %  %*****************************************************************************************************************************
    %       [X_tar, ~]=normalize_points(X_tar);
    %         [X_sub_gt, ~]=normalize_points(X_sub_gt);
    %     X_tar=X_tar';
    %     X_sub_gt=X_sub_gt';
    %           figure(figure)
    %         plot3(X_tar(:,1),X_tar(:,2),X_tar(:,3),'r+','markersize',4); hold on
    % %          view(180,180);
    %         plot3(X_sub_gt(:,1),X_sub_gt(:,2),X_sub_gt(:,3),'bo','markersize',4);
    % %         plot3(vx_affinw(:,1),vx_affinw(:,2),vx_affinw(:,3),'g*','markersize',4);
    %         for i=1:length(X_tar)
    %             plot3([X_sub_gt(i,1);X_tar(i,1)],...
    %                 [X_sub_gt(i,2);X_tar(i,2)],...
    %                 [X_sub_gt(i,3);X_tar(i,3)],'m');
    %         end
    %         xlabel('x'),ylabel('y'),zlabel('z')
    %         axis('equal'); grid on; set (gca, 'box', 'on');
    %         hold off;
    % %    %*****************************************************************************************************************************
    
    
    ind_atlas2validsub_fix=[];
    % do Recog
    output_filename=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\rpm_sc\M\',filename_sub(1:length(filename_sub)-8));
    output_filename=strcat(output_filename,".mat");
    
    [ind_atlas2validsub_pre]=RecogOnPosNoID(X_tar_shape_all,X_sub_shape_all,X_tar, X_sub,arr_pos_var_pwaff,arr_pos_var_pwaff_sc,X_tar_size,X_tar_intensity,X_sub_size,X_sub_intensity,X_sub_gt,output_filename);%noGT
    
    output_filename=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\rpm_sc\rpm_sc\',filename_sub(1:length(filename_sub)-8));
    output_filename=strcat(output_filename,".txt");
    
    fid = fopen(output_filename, 'w');
    for i=1:length(apo_sub_gt)
        if(ind_atlas2validsub_pre(i)<0)
            continue;
        else
%             S=apo_sub{ind_atlas2validsub_pre(i)};
            fprintf(fid, ' %5.3f, %5.3f, %5.3f,%5.3f,0,,,255,255,255\n', X_sub(1,ind_atlas2validsub_pre(i)),...
                X_sub(2,ind_atlas2validsub_pre(i)), X_sub(3,ind_atlas2validsub_pre(i)), index_seg(ind_atlas2validsub_pre(i)));
        end
    end
    fclose(fid);
    output_filename_rename=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\rpm_sc\rpm_sc\',filename_sub(1:length(filename_sub)-8));
    output_filename_rename=strcat(output_filename_rename,".marker");
    copyfile(output_filename, output_filename_rename);
    
    output_filename_swc=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\rpm_sc\rpm_sc_swc\',filename_sub(1:length(filename_sub)-8));
    output_filename=strcat(output_filename_swc,".txt");
    fid = fopen(output_filename, 'w');
    n=0;
    for i=1:length(apo_sub_gt)
        if(ind_atlas2validsub_pre(i)<0)
            continue;
        else
%             S=apo_sub{ind_atlas2validsub_pre(i)};
            n=n+1;
            G=apo_sub_gt{i};
            dis=sqrt((X_sub(1,ind_atlas2validsub_pre(i))-G.x)^2+(X_sub(2,ind_atlas2validsub_pre(i))-G.y)^2+(X_sub(3,ind_atlas2validsub_pre(i))-G.z)^2);
%             if(dis>10)
%                   fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, G.x, G.y, G.z,  n+1);
%             else
            fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, X_sub(1,ind_atlas2validsub_pre(i)),...
                X_sub(2,ind_atlas2validsub_pre(i)), X_sub(3,ind_atlas2validsub_pre(i)),  n+1);
%             end
            n=n+1;
            G=apo_sub_gt{i};
            fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z);
        end
    end
    fclose(fid);
    output_filename_rename=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\rpm_sc\rpm_sc_swc\',filename_sub(1:length(filename_sub)-8));
    output_filename_rename=strcat(output_filename_rename,".swc");
    copyfile(output_filename, output_filename_rename);
     clear X_tar X_sub X_sub_gt
    
end


    function [ind_tar2sub_bi]=RecogOnPosNoID(X_tar_shape_all,X_sub_shape_all,X_tar, X_sub,arr_pos_var_pwaff,arr_pos_var_pwaff_sc,X_tar_size,X_tar_intensity,X_sub_size,X_sub_intensity,X_sub_gt,output_filename)
        ntarcell=length(X_tar);
        nsubcell=length(X_sub);%    min([length(X_sub),length(X_tar)]);
        ind_tar2sub=-1*ones(ntarcell,1);%atlas(i)-->sub(ind_tar2sub(i))
        
        %----------------------------------------------------------------------
        fprintf('\t(1) Do PCA alignment...\n');
        % normaliza point sets
        [X_tar, ~]=normalize_points(X_tar);
        [X_sub, T]=normalize_points(X_sub);
        X_sub_gt=T*[X_sub_gt(1:3,:);ones(1,size(X_sub_gt,2))];
        X_sub_gt=X_sub_gt(1:3,:);
        
        [X_tar_size, ~]=normalize_points(X_tar_size);
        [X_sub_size, ~]=normalize_points(X_sub_size);
        [X_tar_intensity]=normalize_points_1D(X_tar_intensity);
        [X_sub_intensity]=normalize_points_1D(X_sub_intensity);
        for i=1:6
            [X_sub_shape_all{i},~]=normalize_points(X_sub_shape_all{i});
            [X_tar_shape_all{i},~]=normalize_points(X_tar_shape_all{i});
        end
        
        
        % PCA align
        tmp=cov(X_tar'); [~,~,V]=svd(tmp); T_tar=V;
        tmp=cov(X_sub'); [~,~,V]=svd(tmp); T_sub=V;
        T_pca=inv(T_tar')*T_sub'
        X_sub=T_pca*X_sub;
        X_sub_gt=T_pca*X_sub_gt;
        
        %rectify mirror (mirror along 4 dirs and find dir with min dis)
        x=X_sub'; x_bk=x;  yy=X_sub_gt';  x_gt_bk=yy;
        y=X_tar';
        [xmax, dim] = size(x);
        [ymax, dim] = size(y);
        theta = [0,0,0; 0,0,180; 0,180,0; 0,180,180];
        for i=1:4
            sita_x=theta(i,1)/180*pi;  sita_y=theta(i,2)/180*pi;  sita_z=theta(i,3)/180*pi;
            Rx=[1, 0, 0; 0, cos(sita_x), sin(sita_x); 0, -sin(sita_x), cos(sita_x)];
            Ry=[cos(sita_y), 0, -sin(sita_y); 0, 1, 0; sin(sita_y), 0, cos(sita_y)];
            Rz=[cos(sita_z), sin(sita_z), 0; -sin(sita_z), cos(sita_z), 0; 0, 0, 1];
            x=x_bk*Rz*Ry*Rx;
            %记录各情况两点集接近程度tmp(i,j)=dis(y(i),x(j))^2
            tmp = zeros (ymax, xmax);
            for j=1:dim
                tmp = tmp + (y(:,j) * ones(1,xmax) - ones(ymax,1) * x(:,j)').^2;
            end
            near=sort(tmp,2);%sort each rows
            dis(i)=sum(sum(near(:,1:5)));%最近的5个点距离求和
        end
        [q,p]=sort(dis);
        sita_x=theta(p(1),1)/180*pi;  sita_y=theta(p(1),2)/180*pi;  sita_z=theta(p(1),3)/180*pi;
        Rx=[1, 0, 0; 0, cos(sita_x), sin(sita_x); 0, -sin(sita_x), cos(sita_x)];
        Ry=[cos(sita_y), 0, -sin(sita_y); 0, 1, 0; sin(sita_y), 0, cos(sita_y)];
        Rz=[cos(sita_z), sin(sita_z), 0; -sin(sita_z), cos(sita_z), 0; 0, 0, 1];
        x=x_bk*Rz*Ry*Rx;
        X_sub=x';
        yy=x_gt_bk*Rz*Ry*Rx;
        X_sub_gt=yy';
              
        clear i j x y x_bk xmax ymax dim theta i j sita_x sita_y sita_z tmp p q Rx Ry Rz near dis V tmp T_tar T_sub T_pca
        
        %----------------------------------------------------------------------
        fprintf('\t(2) Do RPM matching...\n');
        frac        = 1;
        T_init      = 0.06;%0.006
        T_final     = 0.0001;
        lamda1_init = 0.1;  %0.1 big=affine
        lamda2_init = 0.01; %0.01
        X_sub1=X_sub;
        disp_flag=0;
               
        for iter_rpm=1:1
            
            %             [c,d,vx,m]=cMIX_tps_shape (X_sub',X_tar',frac,T_init, T_final,lamda1_init,lamda2_init,disp_flag,X_sub_gt');
            [c,d,vx,m]=cMIX_tps_shape (X_sub',X_tar',frac,T_init,T_final,lamda1_init,lamda2_init,disp_flag,X_sub_gt',X_tar_shape_all,X_sub_shape_all);
             save (output_filename,'m');
        end
%         load(output_filename,'m');
        P=m;
        for i=1:558
            [maxprob,ind_row]=max(P);%find max in each col
            [~,ind]=max(maxprob);
            ind_tar2sub(ind)=ind_row(ind);
            P(:,ind)=-1; P(ind_row(ind),:)=-1;
        end
        ind_tar2sub_bi=ind_tar2sub;
        clear P
        

    end

    function [xn]=normalize_points_1D(x)
        
        xn=(x-min(x))./(max(x)-min(x));
        
    end
    
    
end


% C.elegant cell recognition based on ARS_segmentation
% by Liyuanyuan 20220309


function main()
clc
clear all
close all
addpath('./RPM/');
addpath('./matlab_io_basicdatatype/');
addpath('./bipartite');

load("data\Statistical_atlas.mat",'average_position', 'APV','RPV','name');

data_input='data\example\input\';
data_gt='data\example\gt\';

L=dir([data_input,'/*.apo']);
for iter_file=1:length(L)
    filename_sub=L(iter_file).name;
    filepath_apo_sub=[data_input,filename_sub];
    filepath_apo_sub_gt=strcat(data_gt,filename_sub(1:length(filename_sub)-8));
    filepath_apo_sub_gt=strcat(filepath_apo_sub_gt,'_mask.apo');
    apo_sub_gt=load_v3d_apo_file(filepath_apo_sub_gt);
        
    fprintf('==========================================================\n');
    fprintf('#(%d):%s\n',iter_file,filename_sub);
    fprintf('==========================================================\n');
    
    fprintf('[1] Load data ...\n');
    % Load sub apo file
    apo_sub=load_v3d_apo_file(filepath_apo_sub);
    fprintf('\tload sub apo file:[%s]\n',filename_sub);
    fprintf('\t\t[%d] cells load\n',length(apo_sub));
    clear ind filepath_apo_sub
    %Load atlas apo file 
    
    fprintf('[4] Do cell recog...\n');
    X_sub_gt=-1*ones(3,558);
    
    for i=1:558
        X_sub_gt(1,i)=apo_sub_gt{i}.x;
        X_sub_gt(2,i)=apo_sub_gt{i}.y;
        X_sub_gt(3,i)=apo_sub_gt{i}.z;    
    end
     
    n=0;
    for i=1:length(apo_sub)  
        n=n+1;
        index_seg(n)=apo_sub{i}.n;
        X_sub(1,n)=apo_sub{i}.x;
        X_sub(2,n)=apo_sub{i}.y;
        X_sub(3,n)=apo_sub{i}.z;     
    end
      
    [ind_atlas2validsub_pre]=RecogOnPosNoID(average_position, X_sub,APV,RPV,X_sub_gt);
    
    output_filename=strcat('data\example\result\',filename_sub(1:length(filename_sub)-8));
    output_filename=strcat(output_filename,".txt");   
    fid = fopen(output_filename, 'w');
     fprintf(fid,"x   ,     y,      z,      pre,   segmentation_label\n");
    for i=1:length(apo_sub_gt)
        if(ind_atlas2validsub_pre(i)<0)
            continue;
        else
            fprintf(fid, ' %5.3f, %5.3f, %5.3f,%5.3f,%s\n',  X_sub(1,ind_atlas2validsub_pre(i)),X_sub(2,ind_atlas2validsub_pre(i)), X_sub(3,ind_atlas2validsub_pre(i)),...
                index_seg(ind_atlas2validsub_pre(i)),name{i});
        end
    end
    fclose(fid);
    
     fprintf('Cell recognition completed !!! \n');
    
%     output_filename_swc=strcat('data\example\result\',filename_sub(1:length(filename_sub)-8));
%     output_filename=strcat(output_filename_swc,".txt");
%     fid = fopen(output_filename, 'w');
%     n=0;
%     for i=1:length(apo_sub_gt)
%         if(ind_atlas2validsub_pre(i)<0)
%             continue;
%         else
%             n=n+1;
%             fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, X_sub(1,ind_atlas2validsub_pre(i)),...
%                 X_sub(2,ind_atlas2validsub_pre(i)), X_sub(3,ind_atlas2validsub_pre(i)),  n+1);
%             n=n+1;
%             G=apo_sub_gt{i};
%             fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z);
%         end
%     end
%     fclose(fid);    
    
    clear X_tar X_sub X_sub_gt
end


    function [ind_tar2sub_bi]=RecogOnPosNoID(X_tar, X_sub,arr_pos_var_pwaff,arr_pos_var_pwaff_sc,X_sub_gt)
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
        disp_flag=0;
           
        [c,d,vx,m]=cMIX_tps_shape (X_sub',X_tar',frac,T_init, T_final,lamda1_init,lamda2_init,disp_flag,X_sub_gt'); 
        
       
        [maxprob,ind_tar2sub]=max(m);
    
        X_sub_or=X_sub;
        X_tar_or=X_tar;
        fprintf('\t(3) Do bipartite cell recog based on APV and RPV ... ...\n');
        for iter=1:25
            ind_tar2sub_valid=find(ind_tar2sub>0);
            X_tar=X_tar_or(:,ind_tar2sub_valid);
            X_sub=X_sub_or(:,ind_tar2sub(ind_tar2sub_valid));%reorder according to current matching result
            X_sub_index=X_sub_or(:,ind_tar2sub(ind_tar2sub_valid));%reorder according to current matching result
            %update ind_tar2sub_fix according to current reorder
            
            % *****************************************************************************************************************************************************************
            % *****************************************************************************************************************************************************************
            if(iter<20)
                T=affine3D_model(X_sub,X_tar);  %T*X_sub=X_tar
                X_sub2tar=T*[X_sub(1:3,:);ones(1,size(X_sub,2))];
                X_sub=X_sub2tar(1:3,:);
                X_sub_or2tar=T*[X_sub_or(1:3,:);ones(1,size(X_sub_or,2))];
                X_sub_all=X_sub_or2tar(1:3,:);
                
                clear T X_sub2tar X_sub_or2tar
                
                %piece-wise affine warp the sub to atlas
                %so that the uneven stretch along worm can be better corrected
          
                %piecewise affine align sub to atlas
                xmin=min(X_sub(1,:));
                xmax=max(X_sub(1,:));
                piecesize=(xmax-xmin)/10;
                piecestep=piecesize/10;
                npiece=0;
                for step=0:100
                    npiece=npiece+1;
                    %find all point within current segment/piece
                    xmin_piece=xmin+step*piecestep;
                    xmax_piece=xmin_piece+piecesize;
                    ind=find(X_sub(1,:)>=xmin_piece & X_sub(1,:)<=xmax_piece);
                    X_sub_piece=X_sub(:,ind);
                    X_tar_piece=X_tar(:,ind);
                    
                    ind_all=find(X_sub_all(1,:)>=xmin_piece & X_sub_all(1,:)<=xmax_piece);
                    X_sub_piece_all=X_sub_all(:,ind_all);
                    
                    if(length(X_tar_piece)<5)
                        continue;
                    end

                    T=affine3D_model(X_sub_piece,X_tar_piece);  %T*X_sub=X_tar
                    X_sub2tar_piece=T*[X_sub_piece(1:3,:);ones(1,size(X_sub_piece,2))];
                    X_sub2tar_piece_all=T*[X_sub_piece_all(1:3,:);ones(1,size(X_sub_piece_all,2))];
              
                    cellarr_X_sub2tar_piece{npiece}=X_sub2tar_piece_all(1:3,:);
                    cellarr_ind_piece{npiece}=ind_all;
                    
                    if(xmax_piece>xmax) break; end
                end
                
                % average point in all pieces to obtain the piecewise affine alignment result
                for i=1:length(X_sub_or)
                    cellarr_avg{i}=[];
                end
                for i=1:length(cellarr_ind_piece)
                    for j=1:length(cellarr_ind_piece{i})
                        ind=cellarr_ind_piece{i}(j);
                        cellarr_avg{ind}=[cellarr_avg{ind},cellarr_X_sub2tar_piece{i}(:,j)];
                    end
                end
                index_valid=[];
                for i=1:length(X_sub_or)
                    value=mean(cellarr_avg{i},2);
                    if( isempty(value))
                        X_sub2tar_avg(:,i)=[0,0,0];
                    else
                        X_sub2tar_avg(:,i)=mean(cellarr_avg{i},2);
                    end
                end
                X_sub2tar_avg_abs=abs(X_sub2tar_avg);
                X_sub2tar_avg_abs=sum(X_sub2tar_avg_abs,1);
                index_value=find(X_sub2tar_avg_abs~=0);
     
                X_sub_index=X_sub_or(:,index_value);
                X_sub=X_sub2tar_avg(:,index_value);
                
                lamda1=0.1;%20*0.5^iter;
                [w,d,k]  = update_tps( X_sub_index', X_sub',lamda1);
                vx = update_marker(X_sub_or',X_sub_index',w,d);
                X_sub_bk_paffine=vx';
                
            else
                lamda1=0.1;%20*0.5^iter;
                [w,d,k]  = update_tps( X_sub_index', X_tar',lamda1);
                vx = update_marker(X_sub_or',X_sub_index',w,d);
                X_sub_bk_paffine=vx';
                
            end
         
            clear i xmin xmax piecesize piecestep npiece step
            clear cellarr_avg cellarr_ind_piece cellarr_X_sub2tar_piece ind j X_sub2tar_avg X_sub2tar_piece
            clear X_tar_piece xmax_piece xmin_piece T X_sub_piece
            
            
            
            % calculate the assignment energy of each atlas point to all topre points
            shape_prob_mse=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
            lam=20;%*0.95^iter;
            
            for i=1:ntarcell
                xyzdiff_atlas2topre=repmat(X_tar_or(:,i),1,nsubcell)-X_sub_bk_paffine;
                exponential_term=0.5*sum(xyzdiff_atlas2topre.^2./(lam^2*repmat(arr_pos_var_pwaff(:,i),1,nsubcell)));
                shape_prob_mse(i,:)=-exponential_term;
            end
            
            [Cost]=ShapeContext3D(X_sub_bk_paffine',X_tar_or');
            
            shape_prob_sc=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
            for i=1:ntarcell
                xyzdiff_atlas2topre=Cost(:,i)';
                exponential_term=0.5*(xyzdiff_atlas2topre.^2./(lam^2*repmat((arr_pos_var_pwaff_sc(:,i)),1,nsubcell)));               
                shape_prob_sc(i,:)=-exponential_term;
            end
        
            clear i  xyzdiff_atlas2topre exponential_term
            
            % do bipartite
            [mat_assignment,assigncost]=munkres(-shape_prob_mse);
            % find atlas to pre matching index
            ind_tar2sub_bi_mse=-1*ones(ntarcell,1);
            for i=1:ntarcell
                assignment=find(mat_assignment(i,:)==1);
                if(~isempty(assignment) )
                    ind_tar2sub_bi_mse(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
                end
            end
            ind_tar2sub_bi_mse_max=-1*ones(ntarcell,1);
            [maxprob,ind_tar2sub_bi_mse_max_col]=max(shape_prob_mse');
            
            [maxprob,ind_tar2sub_bi_mse_max_row]=max(shape_prob_mse);
            for i=1:ntarcell
                if(ind_tar2sub_bi_mse_max_row(ind_tar2sub_bi_mse_max_col(i))==i)
                    ind_tar2sub_bi_mse_max(i)=ind_tar2sub_bi_mse_max_col(i);
                end
            end
           
            [mat_assignment,assigncost]=munkres(-shape_prob_sc);
            % find atlas to pre matching index
            ind_tar2sub_bi_sc=-1*ones(ntarcell,1);
            for i=1:ntarcell
                assignment=find(mat_assignment(i,:)==1);
                if(~isempty(assignment) )
                    ind_tar2sub_bi_sc(i)=assignment;
                end
            end
            ind_tar2sub_bi_sc_max=-1*ones(ntarcell,1);
            [maxprob,ind_tar2sub_bi_sc_max_col]=max(shape_prob_sc');
            
            [maxprob,ind_tar2sub_bi_sc_max_row]=max(shape_prob_sc);
            for i=1:ntarcell
                if(ind_tar2sub_bi_sc_max_row(ind_tar2sub_bi_sc_max_col(i))==i)
                    ind_tar2sub_bi_sc_max(i)=ind_tar2sub_bi_sc_max_col(i);
                end
            end
            
            ind_tar2sub_bi=-1*ones(ntarcell,1);
            for i=1:ntarcell
                n=0;
                select=[ind_tar2sub_bi_mse_max(i),ind_tar2sub_bi_mse(i),ind_tar2sub_bi_sc(i),ind_tar2sub_bi_sc_max(i)];

                [M,F,C] = mode(select);
                if(F>3)
                    ind_tar2sub_bi(i)=M;
                end
            end
            if( length(find(ind_tar2sub_bi ~=-1))<5)
                ind_tar2sub_bi=ind_tar2sub_bi_mse;
                break;
            end
            ind_tar2sub=ind_tar2sub_bi;
        
            clear i mat_assignment assigncost c C Cost d i_shape ind_all ind_tar2sub_bi_mse ind_tar2sub_bi_mse_max ind_tar2sub_bi_mse_max_col ind_tar2sub_bi_mse_max_row
            clear i assignment ind_tar2sub_bi_sc_max  ind_tar2sub_bi_sc_max_col ind_tar2sub_bi_sc_max_row ind_tar2sub_valid index_valid index_value
            clear iter iter_rpm k m M maxprob pca_shape vx w  X_sub1 X_sub2tar_avg_abs X_sub2tar_piece_all  X_sub_index
            clear X_sub_intensity X_sub_or2tar X_sub_piece_all X_sub_shape_all X_tar_intensity X_tar_shape_all yy
            
        end
        
        shape_prob=exp(shape_prob_mse)+exp(shape_prob_sc);
        
        [mat_assignment,assigncost,dMat]=munkres(-shape_prob);
        ind_tar2sub_bi=-1*ones(ntarcell,1);
        for i=1:ntarcell
            assignment=find(mat_assignment(i,:)==1);
            if(~isempty(assignment) )
                ind_tar2sub_bi(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
            end
        end
   
        clear ind_tar2sub ind_tar2sub_bi_sc mat_assignment select shape_prob shape_prob_mse shape_prob_sc
    end
   
end


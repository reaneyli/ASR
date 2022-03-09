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

load('.\traindata\man2atlas_pwaffine_meanstd_norm_sc.mat','arr_pos_mean_pwaff');%arr_pos_mean_pwaff
arr_pos_var_pwaff_sc=arr_pos_mean_pwaff.^2; clear arr_pos_mean_pwaff

%'noID'(output recog result no accu), 'withID'(output recog accu)
datapath='Y:\2.postgrad\16-liyuanyuan\C_elegent\data\wangmeng_pre_mask\segmentation2point_size\';
data_gt="Y:\2.postgrad\16-liyuanyuan\C_elegent\data\cellseganno_4train_20200330\apo_annotation_index_new\"
L=dir([datapath,'/*.apo']);
for iter_file=1:1%length(L)
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
    fprintf('\tload sub apo file:[%s]\n',filename_sub);
    fprintf('\t\t[%d] cells load\n',length(apo_sub));
    clear ind filepath_apo_sub
    %Load atlas apo file
    filepath_apo_atlas='.\traindata\update_atlas_size.apo';
    apo_atlas=load_v3d_pointcloud_file(filepath_apo_atlas);
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
    
    X_sub=-1*ones(3,length(apo_sub));
    X_sub_size=-1*ones(3,length(apo_sub));
    X_sub_intensity=-1*ones(1,length(apo_sub));
    
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
    end
    
    for i=1:length(apo_sub)
        X_sub(1,i)=apo_sub{i}.x;
        X_sub(2,i)=apo_sub{i}.y;
        X_sub(3,i)=apo_sub{i}.z;
        X_sub_size(1,i)=apo_sub{i}.pixmax;
        X_sub_size(2,i)=apo_sub{i}.mass;
        X_sub_size(3,i)=apo_sub{i}.sdev;
        X_sub_intensity(1,i)=apo_sub{i}.volsize;
        
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
    
    [ind_atlas2validsub_pre]=RecogOnPosNoID(X_tar, X_sub,arr_pos_var_pwaff,arr_pos_var_pwaff_sc,X_tar_size,X_tar_intensity,X_sub_size,X_sub_intensity,X_sub_gt);%noGT
    
    output_filename=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\test\rpm_sc\',filename_sub(1:length(filename_sub)-8));
    output_filename=strcat(output_filename,".txt");
    
    fid = fopen(output_filename, 'w');
    for i=1:length(apo_sub_gt)
        S=apo_sub{ind_atlas2validsub_pre(i)};
        fprintf(fid, ' %5.3f, %5.3f, %5.3f,%5.3f,0,,,255,255,255\n', S.x, S.y, S.z,S.n);
    end
    fclose(fid);
    output_filename_rename=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\test\rpm_sc\',filename_sub(1:length(filename_sub)-8));
    output_filename_rename=strcat(output_filename_rename,".marker");
    copyfile(output_filename, output_filename_rename);
    
    output_filename_swc=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\test\rpm_sc_swc\',filename_sub(1:length(filename_sub)-8));
    output_filename=strcat(output_filename_swc,".txt");
    fid = fopen(output_filename, 'w');
    n=0;
    for i=1:length(apo_sub_gt)
        S=apo_sub{ind_atlas2validsub_pre(i)};
        n=n+1;
        fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1);
        n=n+1;
        G=apo_sub_gt{i};
        fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z);
    end
    fclose(fid);
    output_filename_rename=strcat('Y:\2.postgrad\16-liyuanyuan\C_elegent\registration\test\rpm_sc_swc\',filename_sub(1:length(filename_sub)-8));
    output_filename_rename=strcat(output_filename_rename,".swc");
    copyfile(output_filename, output_filename_rename);
    
end


    function [ind_tar2sub_bi]=RecogOnPosNoID(X_tar, X_sub,arr_pos_var_pwaff,arr_pos_var_pwaff_sc,X_tar_size,X_tar_intensity,X_sub_size,X_sub_intensity,X_sub_gt)
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
        
        
        % PCA align
        tmp=cov(X_tar'); [~,~,V]=svd(tmp); T_tar=V;
        tmp=cov(X_sub'); [~,~,V]=svd(tmp); T_sub=V;
        T_pca=inv(T_tar')*T_sub'
        X_sub=T_pca*X_sub;
        X_sub_gt=T_pca*X_sub_gt;
        clear V tmp T_tar T_sub T_pca
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
        
        clear i j x y x_bk xmax ymax dim theta i j sita_x sita_y sita_z tmp p q Rx Ry Rz near dis
        
        %----------------------------------------------------------------------
        fprintf('\t(2) Do RPM matching...\n');
        frac        = 1;
        T_init      = 0.06;%0.006
        T_final     = 0.0001;
        lamda1_init = 0.1;  %0.1 big=affine
        lamda2_init = 0.01; %0.01
        X_sub1=X_sub;
        disp_flag=1;


        for iter_rpm=1:1
            
%             [c,d,vx,m]=cMIX_tps (X_sub',X_sub1',X_tar',frac,T_init,T_final,lamda1_init,lamda2_init,disp_flag,X_tar_size',X_tar_intensity',X_sub_size',X_sub_intensity',X_sub_gt');
%              save m

%             [maxprob,ind_tar2sub]=max(m);
%             ind_tar2sub_valid=find(ind_tar2sub>0);
%             X_tar1=X_tar_prm(:,ind_tar2sub_valid);
%             X_sub1=X_sub_prm(:,ind_tar2sub(ind_tar2sub_valid));%reorder according to current matching result
%             X_sub_index=X_sub1;%reorder according to current matching result
%             %update ind_tar2sub_fix according to current reorder
%             
%             T=affine3D_model(X_sub1,X_tar1);  %T*X_sub=X_tar
%             X_sub2tar1=T*[X_sub1(1:3,:);ones(1,size(X_sub1,2))];
%             X_sub1=X_sub2tar1(1:3,:);
%             clear T X_sub2tar1
%             %piecewise affine align sub to atlas
%             xmin=min(X_sub1(1,:));
%             xmax=max(X_sub1(1,:));
%             piecesize=(xmax-xmin)/10;
%             piecestep=piecesize/10;
%             npiece=0;
%             for step=0:100
%                 npiece=npiece+1;
%                 %find all point within current segment/piece
%                 xmin_piece=xmin+step*piecestep;
%                 xmax_piece=xmin_piece+piecesize;
%                 ind=find(X_sub1(1,:)>=xmin_piece & X_sub1(1,:)<=xmax_piece);
%                 X_sub_piece=X_sub1(:,ind);
%                 X_tar_piece=X_tar1(:,ind);
%                 if(length(X_tar_piece)<5)
%                     continue;
%                 end
%                 
%                 T=affine3D_model(X_sub_piece,X_tar_piece);  %T*X_sub=X_tar
%                 X_sub2tar_piece=T*[X_sub_piece(1:3,:);ones(1,size(X_sub_piece,2))];
%                 
%                 cellarr_X_sub2tar_piece{npiece}=X_sub2tar_piece(1:3,:);
%                 cellarr_ind_piece{npiece}=ind;
%                 
%                 if(xmax_piece>xmax) break; end
%             end
%             
%             % average point in all pieces to obtain the piecewise affine alignment result
%             for i=1:length(X_tar)
%                 cellarr_avg{i}=[];
%             end
%             for i=1:length(cellarr_ind_piece)
%                 for j=1:length(cellarr_ind_piece{i})
%                     ind=cellarr_ind_piece{i}(j);
%                     cellarr_avg{ind}=[cellarr_avg{ind},cellarr_X_sub2tar_piece{i}(:,j)];
%                 end
%             end
%             index_valid=[];
%             for i=1:length(X_tar)
%                 value=mean(cellarr_avg{i},2);
%                 if( isempty(value))
%                     X_sub2tar_avg1(:,i)=0;
%                 else
%                     X_sub2tar_avg1(:,i)=mean(cellarr_avg{i},2);
%                 end
%             end
%             X_sub2tar_avg_abs=abs(X_sub2tar_avg1);
%             X_sub2tar_avg_abs=sum(X_sub2tar_avg_abs,1);
%             index_value=find(X_sub2tar_avg_abs~=0);
%             X_sub=X_sub2tar_avg1(:,index_value);
%             X_sub_index=X_sub_index(:,index_value);
%             lamda1=50;
%             [w,d,k]  = update_tps( X_sub_index', X_sub',lamda1);
%             vx = update_marker(X_sub_prm',X_sub_index',w,d);
%             X_sub1=vx';
            
        end
                load("m.mat",'m');
        
        %         m_or=m;
        %         for i=1:ntarcell
        %             [maxprob,ind_row]=max(m);%find max in each col
        %             [~,ind]=max(maxprob);
        %             ind_tar2sub(ind)=ind_row(ind);
        %             m(:,ind)=-1; m(ind_row(ind),:)=-1;
        %         end
        [maxprob,ind_tar2sub]=max(m);
        
        %           [mat_assignment,assigncost]=munkres(-m');
        %             % find atlas to pre matching index
        %             ind_tar2sub_bi_mse=-1*ones(ntarcell,1);
        %             for i=1:ntarcell
        %                 assignment=find(mat_assignment(i,:)==1);
        %                 if(~isempty(assignment) )
        %                     ind_tar2sub(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
        %                 end
        %             end
        
        X_sub_or=X_sub;
        X_tar_or=X_tar;
        for iter=1:6
            fprintf('\t(3) Affine align sub to atlas according to RPM matching result...\n');
            ind_tar2sub_valid=find(ind_tar2sub>0);
            X_tar=X_tar_or(:,ind_tar2sub_valid);
            X_sub=X_sub_or(:,ind_tar2sub(ind_tar2sub_valid));%reorder according to current matching result
            X_sub_index=X_sub_or(:,ind_tar2sub(ind_tar2sub_valid));%reorder according to current matching result
            %update ind_tar2sub_fix according to current reorder
           if(iter<13) 
            T=affine3D_model(X_sub,X_tar);  %T*X_sub=X_tar
            X_sub2tar=T*[X_sub(1:3,:);ones(1,size(X_sub,2))];
            X_sub=X_sub2tar(1:3,:);
            clear T X_sub2tar
            
            %piece-wise affine warp the sub to atlas
            %so that the uneven stretch along worm can be better corrected
            fprintf('\t(4) Piece-wise affine warp the sub to atlas ...\n');
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
                if(length(X_tar_piece)<5)
                    continue;
                end
                
                %                                 figure(figure);
                %                                 plot3(X_sub_piece(1,:),X_sub_piece(2,:),X_sub_piece(3,:),'r+','markersize',4); hold on
                %                                 plot3(X_tar_piece(1,:),X_tar_piece(2,:),X_tar_piece(3,:),'bo','markersize',4);
                %                                 for i=1:length(X_tar_piece)
                %                                     plot3([X_tar_piece(1,:);X_sub_piece(1,:)],...
                %                                         [X_tar_piece(2,:);X_sub_piece(2,:)],...
                %                                         [X_tar_piece(3,:);X_sub_piece(3,:)],'m');
                %                                 end
                %                                 xlabel('x'),ylabel('y'),zlabel('z')
                %                                 axis('equal'); grid on; set (gca, 'box', 'on');
                
                T=affine3D_model(X_sub_piece,X_tar_piece);  %T*X_sub=X_tar
                X_sub2tar_piece=T*[X_sub_piece(1:3,:);ones(1,size(X_sub_piece,2))];
                
                %                 figure(figure);
                %                 plot3(X_sub_piece(1,:),X_sub_piece(2,:),X_sub_piece(3,:),'r+','markersize',4); hold on
                %                 plot3(X_sub2tar_piece(1,:),X_sub2tar_piece(2,:),X_sub2tar_piece(3,:),'bo','markersize',4);
                %                 xlabel('x'),ylabel('y'),zlabel('z')
                %                 axis('equal'); grid on; set (gca, 'box', 'on');
                
                cellarr_X_sub2tar_piece{npiece}=X_sub2tar_piece(1:3,:);
                cellarr_ind_piece{npiece}=ind;
                
                if(xmax_piece>xmax) break; end
            end

            
            % average point in all pieces to obtain the piecewise affine alignment result
            for i=1:length(X_tar)
                cellarr_avg{i}=[];
            end
            for i=1:length(cellarr_ind_piece)
                for j=1:length(cellarr_ind_piece{i})
                    ind=cellarr_ind_piece{i}(j);
                    cellarr_avg{ind}=[cellarr_avg{ind},cellarr_X_sub2tar_piece{i}(:,j)];
                end
            end
            index_valid=[];
            for i=1:length(X_tar)
                value=mean(cellarr_avg{i},2);
                if( isempty(value))
                    X_sub2tar_avg(:,i)=0;
                else
                    X_sub2tar_avg(:,i)=mean(cellarr_avg{i},2);
                end
            end
            X_sub2tar_avg_abs=abs(X_sub2tar_avg);
            X_sub2tar_avg_abs=sum(X_sub2tar_avg_abs,1);
            index_value=find(X_sub2tar_avg_abs~=0);
            
            X_sub=X_sub2tar_avg(:,index_value);
            X_sub_index=X_sub_index(:,index_value);
            
            lamda1=0.1;
            [w,d,k]  = update_tps( X_sub_index', X_sub',lamda1);
            vx = update_marker(X_sub_or',X_sub_index',w,d);
            X_sub_bk_paffine=vx';
           else
                  lamda1=0.1;%20*0.5^iter;
                [w,d,k]  = update_tps( X_sub_index', X_tar',lamda1);
                vx = update_marker(X_sub_or',X_sub_index',w,d);
                X_sub_bk_paffine=vx';
               
           end
            
            %             figure(figure);
            %             plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'r+','markersize',4); hold on
            %             plot3(X_sub_bk_paffine(1,:),X_sub_bk_paffine(2,:),X_sub_bk_paffine(3,:),'bo','markersize',4);
            %             xlabel('x'),ylabel('y'),zlabel('z')
            %             axis('equal'); grid on; set (gca, 'box', 'on');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear i xmin xmax piecesize piecestep npiece step
            clear cellarr_avg cellarr_ind_piece cellarr_X_sub2tar_piece ind j X_sub2tar_avg X_sub2tar_piece
            clear X_tar_piece xmax_piece xmin_piece T X_sub_piece
            
            fprintf('\t(5) Do bipartite cell recog based on cell relative pos and variation std ...\n');
            
            % calculate the assignment energy of each atlas point to all topre points
            shape_prob_mse=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
            for i=1:ntarcell
                xyzdiff_atlas2topre=repmat(X_tar_or(:,i),1,nsubcell)-X_sub_bk_paffine;
                exponential_term=0.5*sum(xyzdiff_atlas2topre.^2./(20^2*repmat(arr_pos_var_pwaff(:,i),1,nsubcell)));
                shape_prob_mse(i,:)=-exponential_term;
            end
%             shape_prob_mse=exp(shape_prob_mse);
            [Cost]=ShapeContext3D(X_sub_bk_paffine',X_tar_or');
            
            shape_prob_sc=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
            for i=1:ntarcell
                xyzdiff_atlas2topre=Cost(:,i)';
%                 exponential_term=0.5*((xyzdiff_atlas2topre.^2-repmat(arr_pos_var_pwaff_sc(:,i),1,nsubcell)).^2);
                exponential_term=0.5*(xyzdiff_atlas2topre.^2./(20^2*repmat((arr_pos_var_pwaff_sc(:,i)),1,nsubcell)));
                shape_prob_sc(i,:)=-exponential_term;
            end
%             shape_prob_sc=exp(shape_prob_sc);
            
 %**************************************************************************************************        
%  **************************************************************************************************        
%             figure(figure)
%             plot3(X_sub_or(1,:),X_sub_or(2,:),X_sub_or(3,:),'r+','markersize',4); hold on
%             plot3(X_sub_gt(1,:),X_sub_gt(2,:),X_sub_gt(3,:),'bo','markersize',4);
%             for i=1:length(X_sub_gt)
%                 plot3([X_sub_or(1,ind_tar2sub_bi_sc(i));X_sub_gt(1,i)],...
%                     [X_sub_or(2,ind_tar2sub_bi_sc(i));X_sub_gt(2,i)],...
%                     [X_sub_or(3,ind_tar2sub_bi_sc(i));X_sub_gt(3,i)],'m');
%             end
%             xlabel('x'),ylabel('y'),zlabel('z')
%             axis('equal'); grid on; set (gca, 'box', 'on');
%             hold off; drawnow;
%**************************************************************************************************    
%**************************************************************************************************        
         
            
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
            
%   ***********************************************************************************************************     
%   ***********************************************************************************************************                  
%             figure(figure)
%             plot3(X_sub_or(1,:),X_sub_or(2,:),X_sub_or(3,:),'r+','markersize',4); hold on
%             plot3(X_sub_gt(1,:),X_sub_gt(2,:),X_sub_gt(3,:),'bo','markersize',4);
%             for i=1:length(X_sub_gt)
%                 plot3([X_sub_or(1,ind_tar2sub_bi_mse(i));X_sub_gt(1,i)],...
%                     [X_sub_or(2,ind_tar2sub_bi_mse(i));X_sub_gt(2,i)],...
%                     [X_sub_or(3,ind_tar2sub_bi_mse(i));X_sub_gt(3,i)],'m');
%             end
%             xlabel('x'),ylabel('y'),zlabel('z')
%             axis('equal'); grid on; set (gca, 'box', 'on');
%             hold off; drawnow;
%   ***********************************************************************************************************     
%   ***********************************************************************************************************      
            [mat_assignment,assigncost]=munkres(-shape_prob_sc);
            % find atlas to pre matching index
            ind_tar2sub_bi_sc=-1*ones(ntarcell,1);
            for i=1:ntarcell
                assignment=find(mat_assignment(i,:)==1);
                if(~isempty(assignment) )
                    ind_tar2sub_bi_sc(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
                end
            end
%   ***********************************************************************************************************     
%   ***********************************************************************************************************      
%             figure(figure)
%             plot3(X_sub_or(1,:),X_sub_or(2,:),X_sub_or(3,:),'r+','markersize',4); hold on
%             plot3(X_sub_gt(1,:),X_sub_gt(2,:),X_sub_gt(3,:),'bo','markersize',4);
%             for i=1:length(X_sub_gt)
%                 plot3([X_sub_or(1,ind_tar2sub_bi_sc(i));X_sub_gt(1,i)],...
%                     [X_sub_or(2,ind_tar2sub_bi_sc(i));X_sub_gt(2,i)],...
%                     [X_sub_or(3,ind_tar2sub_bi_sc(i));X_sub_gt(3,i)],'m');
%             end
%             xlabel('x'),ylabel('y'),zlabel('z')
%             axis('equal'); grid on; set (gca, 'box', 'on');
%             hold off; drawnow;
%   ***********************************************************************************************************     
%   ***********************************************************************************************************      
            
            ind_tar2sub_bi=-1*ones(ntarcell,1);
            for i=1:ntarcell
                if ind_tar2sub_bi_sc(i)==ind_tar2sub_bi_mse(i)
                    ind_tar2sub_bi(i)=ind_tar2sub_bi_sc(i);
                end
            end
            
            if( length(find(ind_tar2sub_bi ~=-1))<10)
                break;
            end
            ind_tar2sub=ind_tar2sub_bi;
            
    shape_prob=exp(shape_prob_mse+shape_prob_sc);%+exp(shape_prob_mse)+           
            [mat_assignment,assigncost,dMat]=munkres(-shape_prob);
            % find atlas to pre matching index
            ind_tar2sub_bi=-1*ones(ntarcell,1);
            for i=1:ntarcell
                assignment=find(mat_assignment(i,:)==1);
                if(~isempty(assignment) )
                    ind_tar2sub_bi(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
                end
            end

      %   ***********************************************************************************************************     
%   ***********************************************************************************************************      
            figure(figure)
            plot3(X_sub_or(1,:),X_sub_or(2,:),X_sub_or(3,:),'r+','markersize',4); hold on
            plot3(X_sub_gt(1,:),X_sub_gt(2,:),X_sub_gt(3,:),'bo','markersize',4);
            for i=1:length(X_sub_gt)
                plot3([X_sub_or(1,ind_tar2sub_bi(i));X_sub_gt(1,i)],...
                    [X_sub_or(2,ind_tar2sub_bi(i));X_sub_gt(2,i)],...
                    [X_sub_or(3,ind_tar2sub_bi(i));X_sub_gt(3,i)],'m');
            end
            xlabel('x'),ylabel('y'),zlabel('z')
            axis('equal'); grid on; set (gca, 'box', 'on');
            hold off; drawnow;
%   ***********************************************************************************************************     
%   ***********************************************************************************************************            
  %   ***********************************************************************************************************     
%   ***********************************************************************************************************      
 [maxprob,ind_tar2sub_bi]=max(shape_prob');
figure(figure)
plot3(X_sub_or(1,:),X_sub_or(2,:),X_sub_or(3,:),'r+','markersize',4); hold on
plot3(X_sub_gt(1,:),X_sub_gt(2,:),X_sub_gt(3,:),'bo','markersize',4);
for i=1:length(X_sub_gt)
    plot3([X_sub_or(1,ind_tar2sub_bi(i));X_sub_gt(1,i)],...
        [X_sub_or(2,ind_tar2sub_bi(i));X_sub_gt(2,i)],...
        [X_sub_or(3,ind_tar2sub_bi(i));X_sub_gt(3,i)],'m');
end
xlabel('x'),ylabel('y'),zlabel('z')
axis('equal'); grid on; set (gca, 'box', 'on');
hold off; drawnow;
%   ***********************************************************************************************************     
%   ***********************************************************************************************************            
            
            clear i mat_assignment assigncost
            clear i assignment
            
        end
    end

    function [xn]=normalize_points_1D(x)
        
        xn=(x-min(x))./(max(x)-min(x));
        
    end;
    
    
end


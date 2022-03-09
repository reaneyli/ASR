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

    %'noID'(output recog result no accu), 'withID'(output recog accu)
    datapath='.\input_noID\';   str_flagID='noID';
%     datapath='.\input_withID\'; str_flagID='withID';
    
    L=dir([datapath,'/*.txt']);
    recogaccu=-1*ones(length(L),1);
    for iter_file=1:length(L)
        filename_sub=L(iter_file).name;
        filepath_apo_sub=[datapath,filename_sub];
        
        fprintf('==========================================================\n');
        fprintf('#(%d):%s\n',iter_file,filename_sub);
        fprintf('==========================================================\n');
        
        fprintf('[1] Load data ...\n');
        % Load sub apo file
        apo_sub=load_v3d_apo_file(filepath_apo_sub);
        apo_sub_bk=apo_sub;
        fprintf('\tload sub apo file:[%s]\n',filename_sub);
        fprintf('\t\t[%d] cells load\n',length(apo_sub));
        clear ind filepath_apo_sub
        %Load atlas apo file
        filepath_apo_atlas='.\traindata\b_atlas_train100pwaff_20200417.apo';
        apo_atlas=load_v3d_pointcloud_file(filepath_apo_atlas);
        ind=strfind(filepath_apo_atlas,'\');
        filename_atlas=filepath_apo_atlas(ind(end)+1:end);
        fprintf('\tload atlas apo file:[%s]\n',filename_atlas);
        fprintf('\t\t[%d] cells load\n',length(apo_atlas));
        clear ind filepath_apo_atlas filename_atlas
          
        fprintf('[4] Do cell recog...\n');
        nvalidcell=length(apo_sub);
        % reformat valid apo data to arr
        filepath_apo_sub_gt="D:\Vaa3d\vaa3d_tools\released_plugins\v3d_plugins\C_elegent\c_elegent2\main_c_elegent\itera6\C07E3.6_SD1794_20140116_0001\affine_marker_label.txt";
        load('.\traindata\man2atlas_pwaffine_train100meanstd.mat','arr_pos_std_pwaff');
        arr_pos_var_pwaff=arr_pos_std_pwaff; clear arr_pos_std_pwaff
        apo_sub_gt=load_v3d_apo_file(filepath_apo_sub_gt);
        X_sub_gt=-1*ones(3,558);
        X_tar=-1*ones(3,558); X_sub=-1*ones(3,length(apo_sub));
       
        for i=1:558
            X_tar(1,i)=apo_atlas{i}.x;
            X_tar(2,i)=apo_atlas{i}.y;
            X_tar(3,i)=apo_atlas{i}.z;
            X_sub_gt(1,i)=apo_sub_gt{i}.x;
            X_sub_gt(2,i)=apo_sub_gt{i}.y;
            X_sub_gt(3,i)=apo_sub_gt{i}.z;
        end
 
        for i=1:length(apo_sub)
            X_sub(1,i)=apo_sub{i}.x;
            X_sub(2,i)=apo_sub{i}.y;
            X_sub(3,i)=apo_sub{i}.z;
        end
        
 load('.\traindata\man2atlas_affine_meanstd_norm.mat','arr_pos_std_aff');
arr_pos_var_pwaff=arr_pos_std_aff.^2; clear arr_pos_std_aff
arr_pos_var_pwaff_sum=sum(arr_pos_var_pwaff,1);
        
%         load("M_P.mat");
% vx=vx;
% m_or=m;
% m_std=m;
%     % load('tmp.mat');
%     %find matching index
%      ind_tar2sub_std=-1*ones(length(X_tar),1);
%       ind_tar2sub=-1*ones(length(X_tar),1);
%     for i=1:length(X_tar)
%         [maxprob,ind_row]=min(arr_pos_var_pwaff_sum);
%         arr_pos_var_pwaff_sum(ind_row)=100000;
%         [maxprob,ind]=max(m_std(:,ind_row));
%         ind_tar2sub_std(ind_row)=ind;
%         m_std(:,ind_row)=-1; m_std(ind,:)=-1;
%     end
%       for i=1:length(X_tar)
%         [maxprob,ind_row]=max(m_or);%find max in each col
%         [~,ind]=max(maxprob);
%         ind_tar2sub(ind)=ind_row(ind);
%         m_or(:,ind)=-1; m_or(ind_row(ind),:)=-1;
%     end
%      fid = fopen(".\output\result.txt", 'w');
%     n=0;
%     for i=1:length(X_sub_gt)
%         if ind_tar2sub_std(i)==ind_tar2sub(i)
%             S=apo_sub{ind_tar2sub(i)};
%              n=n+1;
%             fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1); 
%             n=n+1;
%             G=apo_sub_gt{i};
%              fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z); 
%         end
%     
%     end
%     fclose(fid);
%     
%      fid_result_tar_pre = fopen(".\output\result_tar_pre.txt", 'w');
%     n=0;
%     for i=1:length(X_tar)
%     S=apo_sub{ind_tar2sub(i)};
%     GG=apo_sub_gt{i};
%     dis=sqrt((S.x-GG.x)^2+(S.y-GG.y)^2+(S.z-GG.z)^2);
%     if dis>5
%          n=n+1;
%         fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1); 
%         n=n+1;
%         G=apo_atlas{i};
%          fprintf(fid, '%d  2  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z); 
%     end
%     
%     end
%     fclose(fid);
%     
%      fid = fopen(".\output\result_tar_gt.txt", 'w');
%     n=0;
%     for i=1:length(X_sub_gt)
%     S=apo_sub_gt{i};
%      n=n+1;
%     fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1); 
%     n=n+1;
%     G=apo_atlas{i};
%      fprintf(fid, '%d  2  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z); 
%     
%     end
%     fclose(fid);
%            
%         save('C.mat','X_tar','X_sub') 
        
        ind_atlas2validsub_fix=[];
  
        % do Recog
        if(strcmpi(str_flagID,'noID'))
            [ind_atlas2validsub_pre,fid_tar2validsub_pre]=RecogOnPosNoID(X_tar, X_sub, ind_atlas2validsub_fix,apo_sub_gt,apo_sub,X_sub_gt);%noGT
        elseif(strcmpi(str_flagID,'withID'))
            [ind_atlas2validsub_pre,fid_tar2validsub_pre]=RecogOnPosNoID(X_tar, X_sub, ind_atlas2validsub_fix,apo_sub_gt,apo_sub,X_sub_gt, ind_atlas2validsub_man);%hasGT
        else
            fprintf('Invalid str_flagID:[%s], return!\n', str_flagID);
        end
        
        
        fprintf('[5] Save recog results...\n');
        %[noGT]print recog result
        apo_sub=apo_sub_bk;
        %reinitialize apo_sub (remove previous name and comment)
        for i=1:length(apo_sub)
            cellname=strtrim(apo_sub{i}.name);
            apo_sub{i}.fidelity=0.0;
            if ~(strcmpi(cellname,'f1') || strcmpi(cellname,'f2') || strcmpi(cellname,'f3') || strcmpi(cellname,'f4') || strcmpi(cellname,'f5') || ...
                    strcmpi(cellname,'f6') || strcmpi(cellname,'f7') || strcmpi(cellname,'f8') || strcmpi(cellname,'f9') || strcmpi(cellname,'f10') || ...
                    contains(cellname,'nouse') || contains(cellname,'NOUSE'))
                apo_sub{i}.name='';
                apo_sub{i}.comment='*AUTOFAIL*';
            end
        end
        %write recong ID and fidelity to apo_sub
        for i=1:length(apo_atlas)
            if(ind_atlas2validsub_pre(i)<1)%cannot find match
                continue;
            end
            ind_apo_sub=ind_validcell(ind_atlas2validsub_pre(i));
            apo_sub{ind_apo_sub}.name=strtrim(apo_atlas{i}.name);
            apo_sub{ind_apo_sub}.comment='*AUTO*';
            apo_sub{ind_apo_sub}.fidelity=fid_tar2validsub_pre(i);
        end
        %restore *FIX* comment to apo_sub
        for i=1:size(ind_atlas2sub_fix,1)
            apo_sub{ind_atlas2sub_fix(i,2)}.comment='*FIX*';
        end
        %save pre result to file 
        output_filename=['.\output\',filename_sub];
        fp = fopen(output_filename, 'w');
        for i=1:length(apo_sub)
            S=apo_sub{i};
            fprintf(fp, '%d,%s,%s,%s,%d,%d,%d,%.2f,%.2f,%.2f,%d,%d,%.4f\n', ...
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
            S.mass, ...
            S.fidelity);           
        end
        fclose(fp);
        fprintf('(%d):%d validcell recog done!\n',iter_file,nvalidcell);
        fprintf('Save predict result to [%s] done.\n', output_filename);
        
       if(strcmpi(str_flagID,'withID'))
            apo_sub=apo_sub_bk;
            %fill matching info matrix (in atlas order) man-vs-pre
            for i=1:558
                info_man_vs_pre{i,1}=strtrim(apo_atlas{i}.name);       %manual name
                if(ind_atlas2validsub_pre(i)<1)
                    info_man_vs_pre{i,2}=-1;                           %apo index
                    info_man_vs_pre{i,3}='';                           %pre name
                    info_man_vs_pre{i,4}='---';                        %---:miss, xxx:wrong
                else    %找到当前atlas name被赋给了哪个cell并填充
                    ind_apo_sub=ind_validcell(ind_atlas2validsub_pre(i));%当前atlas name被赋给的那个cell的index
                    cell_manID=strtrim(apo_sub{ind_apo_sub}.name);%被赋予新name的cell原始manual name
                    for ii=1:558
                        if(strcmpi(cell_manID,strtrim(apo_atlas{ii}.name)))
                            info_man_vs_pre{ii,2}=ind_apo_sub;
                            info_man_vs_pre{ii,3}=strtrim(apo_atlas{i}.name);
                            break;
                        end
                    end
                    if(ind_atlas2validsub_man(i)~=ind_atlas2validsub_pre(i))
                        info_man_vs_pre{i,4}='xxx';
                    end
                end
            end
            %[hasGT]compute recog accu
            nmiss=558-nvalidcell;
            ind_validpre=find(ind_atlas2validsub_man>0);
            ncorrectpre=length(find(ind_atlas2validsub_man(ind_validpre)==ind_atlas2validsub_pre(ind_validpre)));
            nwrongpre=558-nmiss-ncorrectpre;
            recogaccu(iter_file)=ncorrectpre/558;
            fprintf('(%d):%d cell miss, %d/%d cell recog wrong, accuracy=(558-%d-%d)/558=%.4f\n', ...
                iter_file,nmiss,nwrongpre,nvalidcell,nmiss,nwrongpre,recogaccu(iter_file));
            %print recog analysis result to file (compare man and pre)
            %same as atlas cell order
            output_filename=['.\output\',filename_sub,'_recogreport.txt'];
            fp = fopen(output_filename, 'w');
            fprintf(fp,'[%d]:%d cell miss, %d/%d cell recog wrong, accuracy=(558-%d-%d)/558=%.4f\n', ...
                iter_file,nmiss,nwrongpre,nvalidcell,nmiss,nwrongpre,recogaccu(iter_file));
            fprintf(fp,'#no, z, x, y, id_manual, id_auto, errflag\n');
            for i=1:length(apo_atlas)
                if(ind_atlas2validsub_pre(i)<1)%cannot find match
                  fprintf(fp, '%d,  ,  ,  , %s, , ooooooooooooooo->Can not find match!\n', i, strtrim(apo_atlas{i}.name)); 
                  continue;
                end
                if(ind_atlas2validsub_man(i)<1)
                  fprintf(fp, '%d,  ,  ,  , %s, , mmmmmmmmmmmmmmm->This is the real missed cell!\n', i, strtrim(apo_atlas{i}.name)); 
                  continue;
                end
                errflag='';
                if(ind_atlas2validsub_man(i)~=ind_atlas2validsub_pre(i)) errflag='xxxxxxxxxxxxxxx'; end
                fprintf(fp, '%d, %5.3f, %5.3f, %5.3f, %s, %s, %s\n', i, ...
                    apo_sub{info_man_vs_pre{i,2}}.z, apo_sub{info_man_vs_pre{i,2}}.x, apo_sub{info_man_vs_pre{i,2}}.y, ...
                    info_man_vs_pre{i,1}, info_man_vs_pre{i,3}, errflag);
            end
            fclose(fp);
            fprintf('Save analysis result to [%s] done.\n', output_filename);
        end
    end
end


function [ind_tar2sub,fid_tar2sub]=RecogOnPosNoID(X_tar, X_sub, ind_tar2sub_fix ,X_sub_gt,apo_sub,X_sub_gt1,ind_tar2sub_man)
    ntarcell=length(X_tar);
    nsubcell=length(X_sub);%    min([length(X_sub),length(X_tar)]);
    ind_tar2sub=-1*ones(ntarcell,1);%atlas(i)-->sub(ind_tar2sub(i))
    fid_tar2sub=-1*ones(ntarcell,1);
    X_tar_bk=X_tar; X_sub_bk=X_sub;
    
    %whether groundtruth provided (manual anno for accu cpt)
    b_hasGT=1;
    if nargin<7
        b_hasGT=0;
    end
    K=15;
    Nm=findN(X_sub',K);% Nm为（点数，5）的矩阵，第N行表示第N个点最近点的序号；
    Nf=findN(X_tar',K);
    %----------------------------------------------------------------------
    fprintf('\t(1) Do PCA alignment...\n');
    % normaliza point sets
    [X_tar, ~]=normalize_points(X_tar);
    [X_sub, ~]=normalize_points(X_sub);
     
    [X_sub_gt1, ~]=normalize_points(X_sub_gt1);
    % PCA align
    tmp=cov(X_tar'); [~,~,V]=svd(tmp); T_tar=V;
    tmp=cov(X_sub'); [~,~,V]=svd(tmp); T_sub=V;
    T_pca=inv(T_tar')*T_sub'
    X_sub=T_pca*X_sub;
    clear V tmp T_tar T_sub T_pca
    %rectify mirror (mirror along 4 dirs and find dir with min dis)
    x=X_sub'; x_bk=x;
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
    clear i j x y x_bk xmax ymax dim theta i j sita_x sita_y sita_z tmp p q Rx Ry Rz near dis

    %plot cell point
%     figure;
%     plot3(X_tar(1,:),X_tar(2,:),X_tar(3,:),'r+','markersize',3); hold on
%     plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'bo','markersize',3);
%     for i=1:length(X_tar)
%         if(ind_tar2sub_man(i)<=0)
%             plot3(X_tar(1,i),X_tar(2,i),X_tar(3,i),'ms','markersize',6,'markerfacecolor','m');
%         else
%             plot3([X_tar(1,i);X_sub(1,ind_tar2sub_man(i))],...
%                 [X_tar(2,i);X_sub(2,ind_tar2sub_man(i))],...
%                 [X_tar(3,i);X_sub(3,ind_tar2sub_man(i))]);
%         end
%     end
%     xlabel('x'),ylabel('y'),zlabel('z')
%     axis('equal'); grid on; set (gca, 'box', 'on');
%     hold off;

    %----------------------------------------------------------------------
    fprintf('\t(2) Do RPM matching...\n');
    frac        = 1;
    T_init      = 0.006;%0.006
    T_final     = 0.0001;
    lamda1_init = 0.1;  %0.1 big=affine
    lamda2_init = 0.01; %0.01
    disp_flag   = 1;
   
    
     load('.\traindata\man2atlas_pwaffine_meanstd_norm.mat','arr_pos_mean_pwaff');
        arr_pos_var_pwaff=arr_pos_mean_pwaff.^2; clear arr_pos_mean_pwaff
        arr_pos_var_pwaff_sum=sum(arr_pos_var_pwaff,1);
        % calculate the assignment energy of each atlas point to all topre points
%      shape_prob=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
%         for i=1:ntarcell
%             xyzdiff_atlas2topre=repmat(X_tar_bk(:,i),1,nsubcell)-X_sub;
%             exponential_term=0.5*sum(xyzdiff_atlas2topre.^2./(20^2*repmat(arr_pos_var_pwaff(:,i),1,nsubcell)));
%             shape_prob(i,:)=-exponential_term;
%         end
%         shape_prob=exp(shape_prob);
    X_sub1=X_sub;
   
%     [c,d,vx,m]=cMIX_tps (ind_atlas2sub_man,X_sub',X_tar',frac,T_init,T_final,lamda1_init,lamda2_init,disp_flag);  
%  [c,d,vx,m]=cMIX_tps (X_sub1',X_tar',frac,T_init,T_final,lamda1_init,lamda2_init,disp_flag,Nf,Nm,X_sub_gt1',arr_pos_var_pwaff');
%  save m
% %   load("m.mat");    
% %    vx=vx;
%     m=m;
%     m_std=m;
%     for i=1:ntarcell
%         [maxprob,ind_row]=max(m);%find max in each col
%         [~,ind]=max(maxprob);
%         ind_tar2sub(ind)=ind_row(ind);
%         m(:,ind)=-1; m(ind_row(ind),:)=-1;
%     end
%      ind_tar2sub_std=-1*ones(length(X_tar),1);  
%      arr_pos_var_pwaff_sum1=arr_pos_var_pwaff_sum;
%     for i=1:ntarcell
%         [~,ind_row]=min(arr_pos_var_pwaff_sum1);
%         arr_pos_var_pwaff_sum1(ind_row)=100000;
%         [~,ind]=max(m_std(:,ind_row));
%         ind_tar2sub_std(ind_row)=ind;
%         m_std(:,ind_row)=-1; m_std(ind,:)=-1;
%     end
%     X_sub_op = [];
%     X_tar_op = [];
%     number=1;
%     for j=1:length(X_tar)
%          if ind_tar2sub_std(j)==ind_tar2sub(j)
%              X_sub_op(number,1:3) =X_sub(:,ind_tar2sub(j))';
%               X_tar_op(number,1:3) =X_tar(:,j)';      
%               number=number+1;
%          end
%     end
%     
% %     lamda1=0.1;
% %     [w,d,k]  = update_tps( X_sub_op, X_tar_op,lamda1);
% %      vx = update_marker(X_sub',X_sub_op,w,d);
% %      X_sub1=vx';
%        
% %      save('M_P.mat','vx','m') 
%     
%     clear frac T_init T_final lamda1_init lamda2_init disp_flag c d vx
% 
% %     load('tmp.mat');
% %     find matching index
%     
% %         vano中坐标不能有小数
%     fid = fopen(".\output\result.txt", 'w');
%     n=0;
%     for i=1:length(X_sub_gt)
%     S=apo_sub{ind_tar2sub(i)};
%      n=n+1;
%     fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1); 
%     n=n+1;
%     G=X_sub_gt{i};
%      fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z); 
%     
%     end
%     fclose(fid);
%     
%      fid = fopen(".\output\result_tar_gt.txt", 'w');
%     n=0;
%     for i=1:length(X_sub_gt)
%     S=apo_sub{ind_tar2sub(i)};
%      n=n+1;
%     fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1); 
%     n=n+1;
%     G=X_sub_gt{i};
%      fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z); 
%     
%     end
%     fclose(fid);
%     
%      fid = fopen(".\output\result_tar_pre.txt", 'w');
%     n=0;
%     for i=1:length(X_sub_gt)
%     S=apo_sub{ind_tar2sub(i)};
%      n=n+1;
%     fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1); 
%     n=n+1;
%     G=X_sub_gt{i};
%      fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z); 
%     
%     end
%     fclose(fid);
%      
%     fid = fopen(".\output\result_marker.txt", 'w');
%     for i=1:length(X_sub_gt)
%     S=apo_sub{ind_tar2sub(i)};
%     fprintf(fid, ' %5.3f, %5.3f, %5.3f,0,0,,,255,255,255\n', S.x, S.y, S.z);    
%     end
%     fclose(fid);
%     
%     clear i maxprob ind_row ind m 
% 
%     if b_hasGT
%         length(find(ind_tar2sub_man==ind_tar2sub))/ntarcell
%     end
%     
     X_sub_gt_or=-1*ones(3,length(X_sub_gt1));
     X_sub_pre=-1*ones(3,length(apo_sub));
       
    for i=1:length(X_sub_gt1)
        X_sub_gt_or(1,i)=X_sub_gt{i}.x;
        X_sub_gt_or(2,i)=X_sub_gt{i}.y;
        X_sub_gt_or(3,i)=X_sub_gt{i}.z;
    end
    for i=1:length(apo_sub)
        X_sub_pre(1,i)=apo_sub{i}.x;
        X_sub_pre(2,i)=apo_sub{i}.y;
        X_sub_pre(3,i)=apo_sub{i}.z;
    end
   load("C07E3.6_SD1794_20140116_0001_rpm_mse_sc.mat");    
%    vx=vx;
    m=m;
    m_std=m;
    for i=1:ntarcell
        [maxprob,ind_row]=max(m);%find max in each col
        [~,ind]=max(maxprob);
        ind_tar2sub(ind)=ind_row(ind);
        m(:,ind)=-1; m(ind_row(ind),:)=-1;
    end
     fid = fopen(".\output\result_marker.txt", 'w');
        for i=1:length(X_sub_gt1)
            S=apo_sub{ind_tar2sub(i)};
            fprintf(fid, ' %5.3f, %5.3f, %5.3f,0,0,,,255,255,255\n', S.x, S.y, S.z);    
         end
        fclose(fid);
     [X_sub_gt_or_norm, ~]=normalize_points(X_sub_gt_or);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    figure;
    plot3(X_sub_gt_or_norm(1,:),X_sub_gt_or_norm(2,:),X_sub_gt_or_norm(3,:),'r+','markersize',3); hold on
    plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'bo','markersize',3);
    for i=1:length(X_sub_gt1)
        plot3([X_sub_gt_or_norm(1,i);X_sub(1,ind_tar2sub(i))],...
              [X_sub_gt_or_norm(2,i);X_sub(2,ind_tar2sub(i))],...
              [X_sub_gt_or_norm(3,i);X_sub(3,ind_tar2sub(i))],'m');
    end
    xlabel('x'),ylabel('y'),zlabel('z')
    axis('equal'); grid on; set (gca, 'box', 'on');
    hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%     ----------------------------------------------------------------------
    X_sub_or=X_sub;
    X_tar_or=X_tar;
    for iter=1:6
        fprintf('\t(3) Affine align sub to atlas according to RPM matching result...\n');
        ind_tar2sub_valid=find(ind_tar2sub>0);
        X_tar=X_tar_or(:,ind_tar2sub_valid);
        X_sub=X_sub_or(:,ind_tar2sub(ind_tar2sub_valid));%reorder according to current matching result
         X_sub_index=X_sub_or(:,ind_tar2sub(ind_tar2sub_valid));%reorder according to current matching result
        %update ind_tar2sub_fix according to current reorder

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
        piecesize=(xmax-xmin)/8;
        piecestep=piecesize/8;
        npiece=0;
        for step=0:100
            npiece=npiece+1;
            %find all point within current segment/piece
            xmin_piece=xmin+step*piecestep;
            xmax_piece=xmin_piece+piecesize;
            ind=find(X_sub(1,:)>=xmin_piece & X_sub(1,:)<=xmax_piece);
            X_sub_piece=X_sub(:,ind);
            X_tar_piece=X_tar(:,ind);

            %affine align current segment/piece to atlas
            T=affine3D_model(X_sub_piece,X_tar_piece);  %T*X_sub=X_tar
            X_sub2tar_piece=T*[X_sub_piece(1:3,:);ones(1,size(X_sub_piece,2))];

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
        for i=1:length(X_tar)
            X_sub2tar_avg(:,i)=mean(cellarr_avg{i},2);
        end
        X_sub=X_sub2tar_avg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%% piecewise affine    
%     figure;
%     plot3(X_sub_gt_or_norm(1,:),X_sub_gt_or_norm(2,:),X_sub_gt_or_norm(3,:),'r+','markersize',3); hold on
%     plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'bo','markersize',3);
%     for i=1:length(X_sub_gt1)
%         plot3([X_sub_gt_or_norm(1,i);X_sub(1,i)],...
%               [X_sub_gt_or_norm(2,i);X_sub(2,i)],...
%               [X_sub_gt_or_norm(3,i);X_sub(3,i)],'m');
%     end
%     xlabel('x'),ylabel('y'),zlabel('z')
%     axis('equal'); grid on; set (gca, 'box', 'on');
%     hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        
     lamda1=0.1;
    [w,d,k]  = update_tps( X_sub_index', X_sub',lamda1);
     vx = update_marker(X_sub_or',X_sub_index',w,d);
      X_sub_bk_paffine=vx';  
      
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%     figure;
%     plot3(X_sub_bk_paffine(1,:),X_sub_bk_paffine(2,:),X_sub_bk_paffine(3,:),'r+','markersize',3); hold on
%     plot3(X_tar(1,:),X_tar(2,:),X_tar(3,:),'bo','markersize',3);
%       for i=1:length(X_sub_gt1)
%         plot3([X_sub_gt_or_norm(1,i);X_sub_bk_paffine(1,ind_tar2sub(i))],...
%               [X_sub_gt_or_norm(2,i);X_sub_bk_paffine(2,ind_tar2sub(i))],...
%               [X_sub_gt_or_norm(3,i);X_sub_bk_paffine(3,ind_tar2sub(i))],'m');
%     end
%     xlabel('x'),ylabel('y'),zlabel('z')
%     axis('equal'); grid on; set (gca, 'box', 'on');
%     hold off;    
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%       %%% stps after piecewise affine  
%     figure;
%     plot3(X_sub_gt_or_norm(1,:),X_sub_gt_or_norm(2,:),X_sub_gt_or_norm(3,:),'r+','markersize',3); hold on
%     plot3(X_sub_bk_paffine(1,:),X_sub_bk_paffine(2,:),X_sub_bk_paffine(3,:),'bo','markersize',3);
%     for i=1:length(X_sub_gt1)
%         plot3([X_sub_gt_or_norm(1,i);X_sub_bk_paffine(1,ind_tar2sub(i))],...
%               [X_sub_gt_or_norm(2,i);X_sub_bk_paffine(2,ind_tar2sub(i))],...
%               [X_sub_gt_or_norm(3,i);X_sub_bk_paffine(3,ind_tar2sub(i))],'m');
%     end
%     xlabel('x'),ylabel('y'),zlabel('z')
%     axis('equal'); grid on; set (gca, 'box', 'on');
%     hold off;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear i xmin xmax piecesize piecestep npiece step
    clear cellarr_avg cellarr_ind_piece cellarr_X_sub2tar_piece ind j X_sub2tar_avg X_sub2tar_piece
    clear X_tar_piece xmax_piece xmin_piece T X_sub_piece
    


        fprintf('\t(5) Do bipartite cell recog based on cell relative pos and variation std ...\n');
        %load trained cell pos std
        load('.\traindata\man2atlas_pwaffine_meanstd_norm.mat','arr_pos_std_pwaff');
        arr_pos_var_pwaff=arr_pos_std_pwaff.^2; clear arr_pos_std_pwaff

         load('.\traindata\man2atlas_pwaffine_meanstd_norm_sc.mat','arr_pos_std_pwaff',"arr_pos_mean_pwaff");
        arr_pos_var_pwaff_sc=arr_pos_std_pwaff.^2; clear arr_pos_std_pwaff
         arr_pos_var_pwaff_sc_mean=arr_pos_mean_pwaff.^2; clear arr_pos_mean_pwaff


         load('.\traindata\man2atlas_pwaffine_mean_ppos_std_norm2.mat','arr_pos_std_pwaff');
        arr_pos_var_pwaff_ppos=arr_pos_std_pwaff.^2; clear arr_pos_std_pwaff

        shape_prob_ppos=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points

         X_sub_pp=zeros(3,length(X_sub_bk_paffine)); 
          X_tar_pp=zeros(3,558); 
        for j =1:length(X_sub_bk_paffine)
            x=length(X_sub_bk_paffine(X_sub_bk_paffine(1,:)<X_sub_bk_paffine(1,j)));
            X_sub_pp(1,j)=x/(length(X_sub_bk_paffine)-x);
            y=length(X_sub_bk_paffine(X_sub_bk_paffine(2,:)<X_sub_bk_paffine(2,j)));
            X_sub_pp(2,j)=y/(length(X_sub_bk_paffine)-y);
            z=length(X_sub_bk_paffine(X_sub_bk_paffine(3,:)<X_sub_bk_paffine(3,j)));
            X_sub_pp(3,j)=z/(length(X_sub_bk_paffine)-z);
        end 
         for j =1:length(X_tar_or)
             x=length(X_tar_or(X_tar_or(1,:)<X_tar_or(1,j)));
            X_tar_pp(1,j)=x/(558-x);
            y=length(X_tar_or(X_tar_or(2,:)<X_tar_or(2,j)));
            X_tar_pp(2,j)=y/(558-y);
            z=length(X_tar_or(X_tar_or(3,:)<X_tar_or(3,j)));
            X_tar_pp(3,j)=z/(558-z);
        end 

    %     for i=1:ntarcell
    %         xyzdiff_atlas2topre=repmat(X_tar_pp(1,i),1,nsubcell)-X_sub_pp(1,:);
    %         exponential_term=0.5*(xyzdiff_atlas2topre.^2./(20^2*repmat(arr_pos_var_pwaff_ppos(1,i),1,nsubcell)));
    %         shape_prob_ppos(i,:)=-exponential_term;
    %     end
    %       for i=1:ntarcell
    %         xyzdiff_atlas2topre=repmat(X_tar_pp(1,i),1,nsubcell)-X_sub_pp(1,:);
    %         exponential_term=0.5*(xyzdiff_atlas2topre.^2);
    %         shape_prob_ppos(i,:)=-exponential_term;
    %     end
        for i=1:ntarcell
            xyzdiff_atlas2topre=repmat(X_tar_pp(:,i),1,nsubcell)-X_sub_pp(:,:);
            exponential_term=0.5*sum(xyzdiff_atlas2topre.^2./(20^2*repmat(arr_pos_var_pwaff_ppos(:,i),1,nsubcell)));
            shape_prob_ppos(i,:)=-exponential_term;
        end
    %     for i=1:ntarcell
    %         xyzdiff_atlas2topre=repmat(X_tar_pp(:,i),1,nsubcell)-X_sub_pp(:,:);
    %         exponential_term=0.5*sum(xyzdiff_atlas2topre.^2./(20^2*repmat(arr_pos_var_pwaff(:,i),1,nsubcell)));
    %         shape_prob_ppos(i,:)=-exponential_term;
    %     end

            % calculate the assignment energy of each atlas point to all topre points
            shape_prob_mse=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
            for i=1:ntarcell
                xyzdiff_atlas2topre=repmat(X_tar_or(:,i),1,nsubcell)-X_sub_bk_paffine;
                exponential_term=0.5*sum(xyzdiff_atlas2topre.^2./(20^2*repmat(arr_pos_var_pwaff(:,i),1,nsubcell)));
                shape_prob_mse(i,:)=-exponential_term;
            end
    %        shape_prob=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
    %         for i=1:ntarcell
    %             xyzdiff_atlas2topre=repmat(X_tar(:,i),1,nsubcell)-X_sub_bk_paffine;
    %             exponential_term=0.5*sum((xyzdiff_atlas2topre.^2-repmat(arr_pos_var_pwaff(:,i),1,nsubcell)).^2);
    %             shape_prob(i,:)=-exponential_term;
    %         end

            [Cost]=ShapeContext3D(X_sub_bk_paffine',X_tar_or');

%              shape_prob_sc=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
%             for i=1:ntarcell
%                 xyzdiff_atlas2topre=Cost(:,i)';
%                 exponential_term=0.5*(xyzdiff_atlas2topre.^2./(1^20*repmat(sum(arr_pos_var_pwaff_sc(:,i)),1,nsubcell)));
%                 shape_prob_sc(i,:)=-exponential_term;
%             end

    %     shape_prob_sc=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
    %         for i=1:ntarcell
    %             xyzdiff_atlas2topre=Cost(:,i)';
    %             exponential_term=0.5*(((xyzdiff_atlas2topre.^2-repmat(arr_pos_var_pwaff_sc_mean(:,i),1,nsubcell))./(20^2*repmat(sum(arr_pos_var_pwaff_sc_mean(:,i)),1,nsubcell))).^2);
    %             shape_prob_sc(i,:)=-exponential_term;
    %         end
       shape_prob_sc=zeros(ntarcell,nsubcell);%row:assignment energy of one atlas point to all topre points
            for i=1:ntarcell
                xyzdiff_atlas2topre=Cost(:,i)';
                exponential_term=0.5*((xyzdiff_atlas2topre.^2-repmat(arr_pos_var_pwaff_sc(:,i),1,nsubcell)).^2);
                shape_prob_sc(i,:)=-exponential_term;
            end

    %         for i=1:3
    %             xyzdiff_atlas2topre=X_sub(i,:)-X_tar(i,:)';
    %             a=repmat(arr_pos_var_pwaff(i,:),ntarcell,1);
    %             exponential_term=0.5*(xyzdiff_atlas2topre.^2./(20^2*(a)));
    %             shape_prob=shape_prob-exponential_term;
    %         end
    %         
         shape_prob=exp(shape_prob_mse)+exp(shape_prob_sc);%+exp(shape_prob_mse)+
    %      shape_prob_mse=exp(shape_prob_mse);
    %      shape_prob_ppos=exp(shape_prob_ppos);
    %     
    %     sy_shape_prob_mse         = sum (shape_prob_mse') + max(max(shape_prob_mse));
    %     M_shape_prob_mse         = shape_prob_mse ./ (sy_shape_prob_mse');
    %     
    %     sy_shape_prob_ppos         = sum (shape_prob_ppos') + max(max(shape_prob_ppos));
    %     M_shape_prob_ppos         = shape_prob_ppos ./ (sy_shape_prob_ppos');
    %     
    % % %     c=sum(M);
    % % %     shape_prob = M ./  (sum(M)+max(max(M))) ;
    %     shape_prob=M_shape_prob_mse+M_shape_prob_ppos;

            clear i arr_pos_var_pwaff xyzdiff_atlas2topre exponential_term

            % do bipartite
            [mat_assignment,assigncost]=munkres(-shape_prob_mse);
            % find atlas to pre matching index
            ind_tar2sub_bi_mse=-1*ones(ntarcell,1);
            for i=1:ntarcell
                assignment=find(mat_assignment(i,:)==1);
                if(~isempty(assignment) )
                    ind_tar2sub_bi_mse(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
                    fid_tar2sub(i)=shape_prob_sc(i,assignment);
                end
            end
             figure;
            plot3(X_sub_gt_or(1,:),X_sub_gt_or(2,:),X_sub_gt_or(3,:),'r+','markersize',3); hold on
            plot3(X_sub_pre(1,:),X_sub_pre(2,:),X_sub_pre(3,:),'bo','markersize',3);
            for i=1:length(X_sub_gt1)
                    plot3([X_sub_gt_or(1,i);X_sub_pre(1,ind_tar2sub_bi_mse(i))],...
                          [X_sub_gt_or(2,i);X_sub_pre(2,ind_tar2sub_bi_mse(i))],...
                          [X_sub_gt_or(3,i);X_sub_pre(3,ind_tar2sub_bi_mse(i))],'m');
            end
            xlabel('x'),ylabel('y'),zlabel('z')
            axis('equal'); grid on; set (gca, 'box', 'on');
            hold off;
            
            [mat_assignment,assigncost]=munkres(-shape_prob_sc);
            % find atlas to pre matching index
            ind_tar2sub_bi_sc=-1*ones(ntarcell,1);
            for i=1:ntarcell
                assignment=find(mat_assignment(i,:)==1);
                if(~isempty(assignment) )
                    ind_tar2sub_bi_sc(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
                    fid_tar2sub(i)=shape_prob_sc(i,assignment);
                end
            end
            
             figure;
            plot3(X_sub_gt_or(1,:),X_sub_gt_or(2,:),X_sub_gt_or(3,:),'r+','markersize',3); hold on
            plot3(X_sub_pre(1,:),X_sub_pre(2,:),X_sub_pre(3,:),'bo','markersize',3);
            for i=1:length(X_sub_gt1)
                    plot3([X_sub_gt_or(1,i);X_sub_pre(1,ind_tar2sub_bi_sc(i))],...
                          [X_sub_gt_or(2,i);X_sub_pre(2,ind_tar2sub_bi_sc(i))],...
                          [X_sub_gt_or(3,i);X_sub_pre(3,ind_tar2sub_bi_sc(i))],'m');
            end
            xlabel('x'),ylabel('y'),zlabel('z')
            axis('equal'); grid on; set (gca, 'box', 'on');
            hold off;
            
             ind_tar2sub_bi=-1*ones(ntarcell,1);
           for i=1:ntarcell
                if ind_tar2sub_bi_sc(i)==ind_tar2sub_bi_mse(i)
                    ind_tar2sub_bi(i)=ind_tar2sub_bi_sc(i);
                end
           end
             figure;
            plot3(X_sub_gt_or(1,:),X_sub_gt_or(2,:),X_sub_gt_or(3,:),'r+','markersize',3); hold on
            plot3(X_sub_pre(1,:),X_sub_pre(2,:),X_sub_pre(3,:),'bo','markersize',3);
            for i=1:length(X_sub_gt1)
                if ind_tar2sub_bi(i)>0
                    plot3([X_sub_gt_or(1,i);X_sub_pre(1,ind_tar2sub_bi_sc(i))],...
                          [X_sub_gt_or(2,i);X_sub_pre(2,ind_tar2sub_bi_sc(i))],...
                          [X_sub_gt_or(3,i);X_sub_pre(3,ind_tar2sub_bi_sc(i))],'m');
                end
            end
            xlabel('x'),ylabel('y'),zlabel('z')
            axis('equal'); grid on; set (gca, 'box', 'on');
            hold off;
           
                   
            ind_tar2sub=ind_tar2sub_bi;
                 
             [mat_assignment,assigncost,dMat]=munkres(-shape_prob);
            % find atlas to pre matching index
            ind_tar2sub_bi=-1*ones(ntarcell,1);
            for i=1:ntarcell
                assignment=find(mat_assignment(i,:)==1);
                if(~isempty(assignment) )
                    ind_tar2sub_bi(i)=assignment;%ind_tar2sub(ind_tar2sub_valid(assignment));
                    fid_tar2sub(i)=shape_prob(i,assignment);
                end
            end
          

             figure;
            plot3(X_sub_gt_or(1,:),X_sub_gt_or(2,:),X_sub_gt_or(3,:),'r+','markersize',3); hold on
            plot3(X_sub_pre(1,:),X_sub_pre(2,:),X_sub_pre(3,:),'bo','markersize',3);
            for i=1:length(X_sub_gt1)
                    plot3([X_sub_gt_or(1,i);X_sub_pre(1,ind_tar2sub_bi(i))],...
                          [X_sub_gt_or(2,i);X_sub_pre(2,ind_tar2sub_bi(i))],...
                          [X_sub_gt_or(3,i);X_sub_pre(3,ind_tar2sub_bi(i))],'m');
            end
            xlabel('x'),ylabel('y'),zlabel('z')
            axis('equal'); grid on; set (gca, 'box', 'on');
            hold off;
             
            fid = fopen(".\output\result_marker_bi.txt", 'w');
            for i=1:length(X_sub_gt1)
                S=apo_sub{ind_tar2sub_bi(i)};
                fprintf(fid, ' %5.3f, %5.3f, %5.3f,%5.3f,0,,,255,255,255\n', S.x, S.y, S.z,fid_tar2sub(i));    
             end
            fclose(fid);

             fid = fopen(".\output\result.txt", 'w');
        n=0;
        for i=1:length(X_sub_gt)
        S=apo_sub{ind_tar2sub_bi(i)};
         n=n+1;
        fprintf(fid, '%d 3  %5.3f %5.3f %5.3f  2  %d\n',n, S.x, S.y, S.z,  n+1); 
        n=n+1;
        G=X_sub_gt{i};
         fprintf(fid, '%d  3  %5.3f %5.3f  %5.3f 2  -1\n',n, G.x, G.y, G.z); 

        end
        fclose(fid);
 clear i mat_assignment assigncost ind_tar2sub_bi
        if b_hasGT
            length(find(ind_tar2sub_man==ind_tar2sub))/ntarcell
        end

        clear i assignment

    end
    
end


% C.elegant cell recognition based on manual segmentation
% by Lei Qu 20200315
% modified by Lei Qu @ 20210202

function main()
    clc
    clear all
    close all
    addpath('./RPM/');
    addpath('./matlab_io_basicdatatype/');
    addpath('./bipartite');
    
    datapath='.\input\';
    L=dir([datapath,'/*.ano.ano.txt']);
    recogaccu=-1*ones(length(L),1);
    for ii=1:length(L)
        filename_sub=L(ii).name;
        filepath_apo_sub=[datapath,filename_sub];
        
        fprintf('==========================================================\n');
        fprintf('#(%d):%s\n',ii,filename_sub);
        fprintf('==========================================================\n');
        fprintf('[1] Load data ...\n');
        % Load sub apo file
        apo_sub=load_v3d_apo_file(filepath_apo_sub);
        fprintf('\tload sub apo file:[%s]\n',filename_sub);
        fprintf('\t\t[%d] cells load\n',length(apo_sub));
%         clear ind filepath_apo_sub
        %Load atlas apo file
        filepath_apo_atlas='.\traindata\b_atlas_648pwaff_20200302.apo';
        apo_atlas=load_v3d_pointcloud_file(filepath_apo_atlas);
        ind=strfind(filepath_apo_atlas,'\');
        filename_atlas=filepath_apo_atlas(ind(end)+1:end);
        fprintf('\tload atlas apo file:[%s]\n',filename_atlas);
        fprintf('\t\t[%d] cells load\n',length(apo_atlas));
%         clear ind filepath_apo_atlas filename_atlas

        %Assert validcell num is 558 and generate validcell index
        fprintf('[2] Assert valid cell number = 558...\n');
        ind_validcell=[];
        for i=1:length(apo_sub)
            cellname=apo_sub{i}.name;
            if ~(strcmpi(cellname,'f1') || strcmpi(cellname,'f2') || strcmpi(cellname,'f3') || strcmpi(cellname,'f4') || strcmpi(cellname,'f5') || ...
                    strcmpi(cellname,'f6') || strcmpi(cellname,'f7') || strcmpi(cellname,'f8') || strcmpi(cellname,'f9') || strcmpi(cellname,'f10') || ...
                    contains(cellname,'nouse') || contains(cellname,'NOUSE'))
                ind_validcell=[ind_validcell;i];
            end
        end
        if(length(ind_validcell) == 558)
            fprintf('\tAssert success! nvalidcell=558 \n');
        else
            fprintf('\tAssert fail! nvalidcell=%d, should == 558 --> QUIT!\n', nvalidcell);
%             return;
        end
%         clear i cellname

        % generate groundtruth atlas2sub matching index according to manual annotation 
        % Note: only for viusalization, debugging and compute accu
        ind_atlas2sub_man=zeros(558,1);
        for m=1:length(apo_atlas)
            cellname_atlas=strtrim(apo_atlas{m}.name);
            bfind=0;
            for n=1:length(apo_sub)
                cellname=strtrim(apo_sub{n}.name);
                if(strcmpi(cellname_atlas,cellname))
                    ind_atlas2sub_man(m)=n;
                    bfind=1;break;
                end
            end
            if(bfind==0)
                fprintf('\t%s not found in manual apo!\n',cellname_atlas);
%                 return;
            end
        end
%         clear m n nvalidcell bfind cellname cellname_atlas

        fprintf('[3] Do cell recog...\n');
        % reformat valid apo data to arr
        X_tar=zeros(3,558); X_sub=X_tar;
        for i=1:558
            X_tar(1,i)=apo_atlas{i}.x;
            X_tar(2,i)=apo_atlas{i}.y;
            X_tar(3,i)=apo_atlas{i}.z;
            X_sub(1,i)=apo_sub{ind_validcell(i)}.x;
            X_sub(2,i)=apo_sub{ind_validcell(i)}.y;
            X_sub(3,i)=apo_sub{ind_validcell(i)}.z;
        end
        % update ind_atlas2sub_man to valid cell (only for viusalization, debugging and compute accu)
        ind_atlas2validsub_man=zeros(558,1);
        for i=1:558
            ind_atlas2validsub_man(i)=find(ind_validcell==ind_atlas2sub_man(i));
        end
%         %plot cell point
%         figure;
%         plot3(X_tar(1,:),X_tar(2,:),X_tar(3,:),'r+','markersize',3); hold on
%         plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'bo','markersize',3);
%         for i=1:length(X_tar)
%             plot3([X_tar(1,i);X_sub(1,ind_atlas2validsub_man(i))],...
%                 [X_tar(2,i);X_sub(2,ind_atlas2validsub_man(i))],...
%                 [X_tar(3,i);X_sub(3,ind_atlas2validsub_man(i))]);
%         end
%         xlabel('x'),ylabel('y'),zlabel('z')
%         axis('equal'); grid on; set (gca, 'box', 'on');
%         hold off;
        %do Recog
        ind_atlas2validsub_pre=RecogOnPosNoID(X_tar, X_sub, ind_atlas2validsub_man);
        
        %compute recog accu
        nwrong=length(find(ind_atlas2validsub_man~=ind_atlas2validsub_pre));
        recogaccu(ii)=1-nwrong/length(apo_atlas);
        recogaccu_filename{ii}=filename_sub;
        
        %print recog analysis result to file (compare man and pre)
        %same as atlas cell order
        output_filename=['.\output\',filename_sub(1:end-length('.ano.ano.txt')),'_recog.txt'];
        fp = fopen(output_filename, 'w');
        fprintf(fp,'[%d]:%d/558 cell recog wrong, accuracy=%.4f\n',ii,nwrong,recogaccu(ii));
        fprintf(fp,'#no, z, x, y, id_manual, id_auto, errflag\n');
        for i=1:length(apo_atlas)
            apo_sub_ind=ind_validcell(ind_atlas2validsub_man(i));
            apo_sub_id_man=apo_atlas{i}.name;
            apo_sub_id_pre=apo_atlas{ind_atlas2validsub_man==ind_atlas2validsub_pre(i)}.name;
            errflag='';
            if(ind_atlas2validsub_man(i)~=ind_atlas2validsub_pre(i)) errflag='xxxxxxxxxxxxxxx'; end
            fprintf(fp, '%d, %5.3f, %5.3f, %5.3f, %s, %s, %s\n', i, ...
                apo_sub{apo_sub_ind}.z, apo_sub{apo_sub_ind}.x, apo_sub{apo_sub_ind}.y, ...
                strtrim(apo_sub_id_man), strtrim(apo_sub_id_pre), errflag);
        end
        fclose(fp);
        fprintf('\t(%d):%d/558 cell recog wrong, accuracy=%.4f\n',ii,nwrong,recogaccu(ii));
        fprintf('Save analysis result to [%s] done.\n', output_filename);
    end
end


function ind_tar2sub=RecogOnPosNoID(X_tar, X_sub, ind_tar2sub_man)
    ncell=length(X_tar);
    ind_tar2sub=zeros(ncell,1);%atlas(i)-->sub(ind_tar2sub(i))
    X_tar_bk=X_tar; X_sub_bk=X_sub;
    
    %----------------------------------------------------------------------
    fprintf('\t(1) Do PCA alignment...\n');
    % normaliza point sets
    [X_tar, ~]=normalize_points(X_tar);
    [X_sub, ~]=normalize_points(X_sub);
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
        tmp = zeros (xmax, ymax);
        for j=1:dim
            tmp = tmp + (y(:,j) * ones(1,ymax) - ones(xmax,1) * x(:,j)').^2;
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

%     %plot cell point
%     figure;
%     plot3(X_tar(1,:),X_tar(2,:),X_tar(3,:),'r+','markersize',3); hold on
%     plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'bo','markersize',3);
%     for i=1:length(X_tar)
%         plot3([X_tar(1,i);X_sub(1,ind_tar2sub_man(i))],...
%             [X_tar(2,i);X_sub(2,ind_tar2sub_man(i))],...
%             [X_tar(3,i);X_sub(3,ind_tar2sub_man(i))]);
%     end
%     xlabel('x'),ylabel('y'),zlabel('z')
%     axis('equal'); grid on; set (gca, 'box', 'on');
%     hold off;

    %----------------------------------------------------------------------
    fprintf('\t(2) Do RPM matching...\n');
    frac        = 1;
    T_init      = 0.006;%0.006
    T_final     = 0.0005;
    lamda1_init = 0.1;  %0.1 big=affine
    lamda2_init = 0.01; %0.01
    disp_flag   = 0;
%     [c,d,vx,m]=cMIX_tps (ind_atlas2sub_man,X_sub',X_tar',frac,T_init,T_final,lamda1_init,lamda2_init,disp_flag);
    [c,d,vx,m]=cMIX_tps (X_sub',X_tar',frac,T_init,T_final,lamda1_init,lamda2_init,disp_flag);
    clear frac T_init T_final lamda1_init lamda2_init disp_flag c d vx

    % load('tmp.mat');
    %find matching index
    for i=1:ncell
        [maxprob,ind_row]=max(m);
        [~,ind]=max(maxprob);
        ind_tar2sub(ind)=ind_row(ind);
        m(:,ind)=-1; m(ind_row(ind),:)=-1;
    end
    clear i maxprob ind_row ind m 

    1-length(find(ind_tar2sub_man~=ind_tar2sub))/ncell

    %----------------------------------------------------------------------
    % load('tmp.mat');
    for iter=1:3
        fprintf('\t(3) Affine align sub to atlas according to RPM matching result...\n');
        X_tar=X_tar_bk;
        X_sub=X_sub_bk(:,ind_tar2sub);

        T=affine3D_model(X_sub,X_tar);  %T*X_sub=X_tar
        X_sub2tar=T*[X_sub(1:3,:);ones(1,size(X_sub,2))];
        X_sub=X_sub2tar(1:3,:);
        clear T X_sub2tar

%         figure;
%         plot3(X_tar(1,:),X_tar(2,:),X_tar(3,:),'r+','markersize',3); hold on
%         plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'bo','markersize',3);
%         for i=1:length(X_tar)
%             plot3([X_tar(1,i);X_sub(1,i)],...
%                   [X_tar(2,i);X_sub(2,i)],...
%                   [X_tar(3,i);X_sub(3,i)]);
%         end
%         xlabel('x'),ylabel('y'),zlabel('z')
%         axis('equal'); grid on; set (gca, 'box', 'on');
%         hold off;

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
        for i=1:ncell
            cellarr_avg{i}=[];
        end
        for i=1:length(cellarr_ind_piece)
            for j=1:length(cellarr_ind_piece{i})
                ind=cellarr_ind_piece{i}(j);
                cellarr_avg{ind}=[cellarr_avg{ind},cellarr_X_sub2tar_piece{i}(:,j)];
            end
        end
        for i=1:ncell
            X_sub2tar_avg(:,i)=mean(cellarr_avg{i},2);
        end
        X_sub=X_sub2tar_avg;

%         figure;
%         plot3(X_tar(1,:),X_tar(2,:),X_tar(3,:),'r+','markersize',3); hold on
%         plot3(X_sub(1,:),X_sub(2,:),X_sub(3,:),'bo','markersize',3);
%         for i=1:length(X_tar)
%             plot3([X_tar(1,i);X_sub(1,i)],...
%                   [X_tar(2,i);X_sub(2,i)],...
%                   [X_tar(3,i);X_sub(3,i)]);
%         end
%         xlabel('x'),ylabel('y'),zlabel('z')
%         axis('equal'); grid on; set (gca, 'box', 'on');
%         hold off;

        clear i xmin xmax piecesize piecestep npiece step
        clear cellarr_avg cellarr_ind_piece cellarr_X_sub2tar_piece ind j X_sub2tar_avg X_sub2tar_piece
        clear X_tar_piece xmax_piece xmin_piece T X_sub_piece

        fprintf('\t(5) Do bipartite cell recog based on cell relative pos and variation std ...\n');
        %load trained cell pos std
        load('.\traindata\man2atlas_pwaffine_meanstd.mat','arr_pos_std_pwaff');
        arr_pos_var_pwaff=arr_pos_std_pwaff.^2; clear arr_pos_std_pwaff
        % calculate the assignment energy of each atlas point to all topre points
        shape_prob=zeros(ncell,ncell);%row:assignment energy of one atlas point to all topre points
        for i=1:ncell
            xyzdiff_atlas2topre=repmat(X_tar(:,i),1,ncell)-X_sub;
            exponential_term=0.5*sum(xyzdiff_atlas2topre.^2./(20^2*repmat(arr_pos_var_pwaff(:,i),1,ncell)));
            shape_prob(i,:)=-exponential_term;
        end
        shape_prob=exp(shape_prob);
        clear i arr_pos_var_pwaff xyzdiff_atlas2topre exponential_term
        % do bipartite
        [mat_assignment,assigncost]=munkres(-shape_prob);
        % find atlas to pre matching index
        ind_tar2sub_bi=ind_tar2sub;
        for i=1:ncell
            assignment(i)=find(mat_assignment(i,:)==1);
            ind_tar2sub_bi(i)=ind_tar2sub(assignment(i));
        end
        ind_tar2sub=ind_tar2sub_bi;
        clear i mat_assignment assigncost ind_tar2sub_bi

        1-length(find(ind_tar2sub_man~=ind_tar2sub))/ncell

        clear i assignment

    end
    
end


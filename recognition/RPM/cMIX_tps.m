function [o1,o2,o3,o4] =  cMIX_tps (X_sub,in1,in2,in3,in4,in5,in6,in7,in8,X_tar_size,X_tar_intensity,X_sub_size,X_sub_intensity,X_sub_gt)
% function [o1,o2,o3,o4] =  cMIX_tps (apo_topre_reorder_ind,in1,in2,in3,in4,in5,in6,in7,in8)

if in8
    h=figure;
end
% fig(1); clf; whitebg('k'); set(gcf,'color',[0 0 0]);
% set(gcf,'DoubleBuffer','on')

% Init control parameters:
perT_maxit  = 10;
anneal_rate = 0.95;
SC_T_init      = 0.05;%0.006
anneal_rate_SC = 1.001;
T_STPS=200;
T_std= 0.08;
anneal_rate_std=0.95;

lamda1_init = in6;
lamda2_init = in7;
disp_flag = in8;
x       = in1;
y       = in2;
frac	= in3;
T_init	= in4;
T_final = in5;
trans_type = 'tps';
sigma     = 1;
z         = x;

% init x,y,z:
[xmax, dim] = size(x); x = x (1:frac:xmax, :); [xmax, dim] = size(x);
[ymax, tmp] = size(y); y = y (1:frac:ymax, :); [ymax, tmp] = size(y);
[zmax, tmp] = size(z);
z = x;

% init m:
m              = ones (xmax, ymax) ./ (xmax * ymax);

T0             = max(x(:,1))^2;
moutlier       = 1/sqrt(T0)*exp(-1);       % /xmax *0.001;
% moutlier       = 1/xmax*0.01;
m_outliers_row = ones (1,ymax) * moutlier;
m_outliers_col = ones (xmax,1) * moutlier;

m_outliers_row_intenstity=ones (1,ymax) *1/sqrt(max(X_sub_intensity(:,1))^2)*exp(-1);

m_outliers_row_size=ones (1,ymax) *1/sqrt(max(X_sub_size(:,1))^2)*exp(-1);

% init transformation parameters:
theta = 0; t = zeros (2,1); s = 1;
c_tps = zeros (xmax,dim+1);
d_tps = eye   (dim+1, dim+1);
w     = zeros (xmax+dim+1, dim);

% -------------------------------------------------------------------
% Annealing procedure:
% -------------------------------------------------------------------
T       = T_init;
T_SC      = SC_T_init;
lamda1_stps=50;

vx = x;
vy = y;
vx_affinw=x;

it_total = 1;
flag_stop = 0;
% if disp_flag
%     figure(h)
%     plot3(X_sub_gt(:,1),X_sub_gt(:,2),X_sub_gt(:,3),'r+','markersize',4); hold on
%     plot3(X_sub(:,1),X_sub(:,2),X_sub(:,3),'bo','markersize',4);
%     xlabel('x'),ylabel('y'),zlabel('z')
%     axis('equal'); grid on; set (gca, 'box', 'on');
% end

% options=struct('r_max',2,'r_min',1e-3,'r_bins',15,'a_bins',15,'rotate',0,'method',1,'maxdist',5);
% Normals2=[];
% F2=getHistogramFeatures(y,y,Normals2,options);

while (flag_stop ~= 1)
    for i=1:perT_maxit     % repeat at each termperature.
        % Given vx, y, Update m:
        m = cMIX_calc_m (vx, y, T, m_outliers_row, m_outliers_col,T_SC,X_tar_size,X_tar_intensity,X_sub_size,X_sub_intensity,m_outliers_row_intenstity,m_outliers_row_size);
        
        % Given m, update transformation:
        %         a= (sum(m'))' * ones(1,558);
        %         vy = m * y ./ ( (sum(m'))' * ones(1,dim));
        %
        %         vy1= m  ./ ( (sum(m'))' * ones(1,558))* y;
        vy= m * y;
        
        lamda1 = lamda1_init*length(x)*T;
        lamda2 = lamda2_init*length(x)*T;
        
        
        [c_tps, d_tps, w] = cMIX_calc_transformation (trans_type, ...
            lamda1, lamda2, sigma, x, vy, z);
        [vx_affinw,vx] = cMIX_warp_pts (trans_type, x, z, c_tps, d_tps, w, sigma);
        
        %             vx=pieceaffine(vx',vy');
%          lamda1 = lamda1_stps*length(x)*T;
%                 [w,d,k]  = update_tps(x, vy, lamda1);
%                 % Update the warping template x^w
%                 vx = update_xw(x,w,d,k);
        
        
    end  % end of iteration/perT
    
    T = T * anneal_rate;
    T_SC=T_SC*anneal_rate_SC;
    T_STPS=T_STPS*0.98;
    T_std=T_std*anneal_rate_std;
    
    % Determine if it's time to stop:
    %     fprintf ('T = %.4f:\t lamda1: %.4f lamda2: %.4f\n', T, lamda1, lamda2);
    if T < T_final; flag_stop = 1; end;
    
    % Display:
    if disp_flag
        P=m;
%         for i=1:558
%             [maxprob,ind_row]=max(P);%find max in each col
%             [~,ind]=max(maxprob);
%             ind_tar2sub(ind)=ind_row(ind);
%             P(:,ind)=-1; P(ind_row(ind),:)=-1;
%         end
                   [maxprob,ind_tar2sub]=max(P);
        
        
        it_total = it_total + 1;
        figure(h)
        plot3(X_sub_gt(:,1),X_sub_gt(:,2),X_sub_gt(:,3),'r+','markersize',4); hold on
        %          view(180,180);
        plot3(X_sub(:,1),X_sub(:,2),X_sub(:,3),'bo','markersize',4);
        plot3(y(:,1),y(:,2),y(:,3),'g*','markersize',4);
        for i=1:length(y)
            plot3([X_sub(ind_tar2sub(i),1);X_sub_gt(i,1)],...
                [X_sub(ind_tar2sub(i),2);X_sub_gt(i,2)],...
                [X_sub(ind_tar2sub(i),3);X_sub_gt(i,3)],'m');
        end
        
         for i=1:length(y)
            plot3([X_sub_gt(i,1);y(i,1)],...
                [X_sub_gt(i,2);y(i,2)],...
                [X_sub_gt(i,3);y(i,3)],'g');
        end
        
        xlabel('x'),ylabel('y'),zlabel('z')
        axis('equal'); grid on; set (gca, 'box', 'on');
              
        hold off; drawnow;
    end
end % end of annealing.

o1 = c_tps;
o2 = d_tps;
o3 = vx;
o4 = m;
end

function [m, m_outliers_row, m_outliers_col] = cMIX_calc_m ...
    (vx, y, T, m_outliers_row, m_outliers_col,T_local,X_tar_size,X_tar_intensity,X_sub_size,X_sub_intensity,m_outliers_row_intenstity,m_outliers_row_size)

[xmax,dim] = size(vx);
[ymax,dim] = size(y);

% [Cost]=ShapeContext3D(vx,y);
%   [Cost_std]=StdCost(vx_affinw,y,arr_pos_var_pwaff);

% Given v=tranformed(x), update m:
y_tmp = zeros (xmax, ymax);

for it_dim=1:dim
    y_tmp=y_tmp+ (vx(:,it_dim) - y(:,it_dim)').^2;%-(T_std*ones(xmax,1) * arr_pos_var_pwaff(:,it_dim)' ));
end

% y_size = zeros (xmax, ymax);
% for it_dim=1:dim
%     y_size=y_size+ (X_sub_size(:,it_dim) - X_tar_size(:,it_dim)').^2;
% end
% 
% y_intensity = zeros (xmax, ymax);
% for it_dim=1:1
%     y_intensity=y_intensity+ (X_sub_intensity(:,it_dim) - X_tar_intensity(:,it_dim)').^2;
% end


% m_coor = 1/sqrt(T*T_local) .* exp (-y_tmp/T-Cost/T_local);
m_coor = 1/sqrt(T*T) .* exp (-y_tmp/T);
% m_coor = 1/sqrt(T) .* exp (-y_tmp/T);
% normalize accross the outliers as well:
sy_coor         = sum (m_coor) + m_outliers_row;
m_coor         = m_coor ./ (ones(xmax,1) * sy_coor);
m_coor         =m_coor  ./ ( (sum(m_coor'))' * ones(1,ymax));

%  m_size = 1/sqrt(T) .* exp (-y_size/T);
% sy_size        = sum (m_size) + m_outliers_row_size;
% m_size         = m_size ./ (ones(xmax,1) * sy_size);
% m_size         =m_size  ./ ( (sum(m_size'))' * ones(1,ymax));
% 
%  m_intensity = 1/sqrt(T) .* exp (-y_intensity/T);
% sy_intensity        = sum (m_intensity) + m_outliers_row_intenstity;
% m_intensity         = m_intensity ./ (ones(xmax,1) * sy_intensity);
% m_intensity         =m_intensity  ./ ( (sum(m_intensity'))' * ones(1,ymax));
% 
% M=0.6*m_coor+0.1*m_size+0.3*m_intensity;
% sy_M    = sum (M) ;
% m_M        = M ./ (ones(xmax,1) * sy_M);
% m_M         =m_M  ./ ( (sum(m_M'))' * ones(1,ymax));

%  m_coor = 1/sqrt(T*T_local) .* exp (-y_tmp/T-y_intensity/T_local);
% sy_coor         = sum (m_coor) + m_outliers_row;
% m_coor         = m_coor ./ (ones(xmax,1) * sy_coor);
% m_coor         =m_coor  ./ ( (sum(m_coor'))' * ones(1,ymax));


m=m_coor;
end


function [vx]=pieceaffine(X_sub,X_tar)
[dim,nsubcell] = size(X_sub);
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
for i=1:nsubcell
    cellarr_avg{i}=[];
end
for i=1:length(cellarr_ind_piece)
    for j=1:length(cellarr_ind_piece{i})
        ind=cellarr_ind_piece{i}(j);
        cellarr_avg{ind}=[cellarr_avg{ind},cellarr_X_sub2tar_piece{i}(:,j)];
    end
end
for i=1:nsubcell
    X_sub2tar_avg(:,i)=mean(cellarr_avg{i},2);
end
vx=X_sub2tar_avg';
end

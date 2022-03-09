function F=getHistogramFeatures(Points,Pointsl,Normals,options)
rlog_max=log(options.r_max);
rlog_min=log(options.r_min);
switch size(Points,2)
    case 2
        F=zeros(options.a_bins*options.r_bins,size(Points,1));
        for i=1:size(Points,1)
            % The Current Point
            P=Points(i,:);
            
            % Determine the log-distance an angle of all points
            % relative to the current point
            O=bsxfun(@minus,Pointsl,P);
            if(~isempty(Normals))
                N=Normals(i,:);
                M=[N(:),[-N(2);N(1)]];
                O=O*M;
%                 figure, hold on;
%                 plot(O(:,1),O(:,2),'b.'), 
%                 plot(O(i,1),O(i,2),'r*');
%                 N=N*M;
%                 plot(O(i,1)+N(1,1),O(i,2)+N(1,2),'r*');
%                 pause;
            end
            
            R=sqrt(sum(O.^2,2));
            A = (atan2(O(:,1),O(:,2))+pi)/(2*pi);
            R(R<options.r_min)=options.r_min;
            Rlog = log(R);
            
            % Scale radius and angle between 0..1
            Rlog=(Rlog-rlog_min)/(rlog_max-rlog_min);
            
            % Histogram positions
            Abin=floor(A*options.a_bins);
            Rbin=floor(Rlog*options.r_bins);
            
            % Postions must be inside the histogram
            Abin(Abin>options.a_bins-1)=options.a_bins-1; Abin(Abin<0)=0;
            Rbin(Rbin>options.r_bins-1)=options.r_bins-1; Rbin(Rbin<0)=0;
            
            % Construct the polar-distance histogram by counting
            %x=(1:size(Pointsl,1))';
            y=Abin*options.r_bins+Rbin+1;
            %ind=x+(y-1).*size(Pointsl,1);
            %H=zeros(size(Pointsl,1),options.a_bins*options.r_bins);
            %H(ind)=1;
            %H=sum(H,1)./size(Pointsl,1);
            H=accumarray(y,1,[options.a_bins*options.r_bins 1]);
            H=H./size(Pointsl,1);
            % Reshape the histogram vector to the a 2D histogram image
            % H=reshape(H,[options.r_bins,options.a_bins]);
            
            % Store the Histogram
            F(:,i)=H(:);
        end
    case 3
        options.w_bins=options.a_bins;
        F=zeros(options.w_bins*options.a_bins*options.r_bins,size(Points,1));
%         Nm=findN(Points,100);% Nm为（点数，5）的矩阵，第N行表示第N个点最近点的序号；
        for i=1:size(Points,1)
%             Pointsl=Points(Nm(i,:),:);
            % The Current Point
            P=Points(i,:);     
            % Determine the log-distance an angle of all points
            % relative to the current point
            O=Pointsl-P;%bsxfun(@minus,Pointsl,P);
            R=sqrt(sum(O.^2,2));
            A = (atan2(O(:,1),O(:,2))+pi)/(2*pi);
            R(R<options.r_min)=options.r_min;
            W = ((O(:,3)./R)+1)/2;
            Rlog = log(R);
            
            % Scale radius and angle between 0..1
            Rlog=(Rlog-rlog_min)/(rlog_max-rlog_min);
            
            % Histogram positions
            Abin=A*(options.a_bins-1e-10);
            Wbin=W*(options.w_bins-1e-10);
            Rbin=Rlog*(options.r_bins-1e-10);
            
			% Construct the polar-distance histogram by counting
			Abin=floor(Abin); Wbin=floor(Wbin); Rbin=floor(Rbin);
			y= Wbin*options.r_bins*options.a_bins + Abin*options.r_bins + Rbin + 1;
            % This is slow ....
            %x=(1:size(Pointsl,1))';
            %ind=x+(y-1).*size(Pointsl,1);
            %H=false(size(Pointsl,1),options.w_bins*options.a_bins*options.r_bins);
			%H(ind)=true;
            %H=sparse(x,y,1,size(Pointsl,1),options.w_bins*options.a_bins*options.r_bins,size(Pointsl,1));
            %H=sum(H,1)./size(Pointsl,1);
            H=accumarray(y,1,[options.w_bins*options.a_bins*options.r_bins 1]);
            H=H./size(Pointsl,1);
            % Reshape the histogram vector to the a 2D histogram image
            %H=reshape(H,[options.r_bins,options.a_bins, options.w_bins]);
            
            % Store the Histogram
            F(:,i)=H(:);
        end
end

function [Pointsn,Pointsln]=NormalizePoints(Points,Pointsl,options)
Pointst=bsxfun(@minus,Points,mean(Points,1));
Pointslt=bsxfun(@minus,Pointsl,mean(Points,1));
if(size(Points,2)==2) % 2D or 3D
    if(options.rotate>1)
        % Eigen Vector align
        [Evalues, Evectors]=PCA(Pointst');
        
        % Do the Point transformation
        M=eye(2)/Evectors;
        Ml=zeros(3,3); Ml(1:2,1:2)=M;  Ml(3,3)=1;
        Pointsr=transformvertices(Pointst,Ml);
		Pointslr=transformvertices(Pointslt,Ml);
		
        % Heaviest axis flip
        R=sqrt(sum(Pointsr.^2,2));
        h=Pointsr>0;
        m(1)=sum(R(h(:,1)));
        m(2)=sum(R(~h(:,1)));
        m(3)=sum(R(h(:,2)));
        m(4)=sum(R(~h(:,2)));
        if(m(2)>m(1)), Pointsr(:,1)=-Pointsr(:,1); Pointslr(:,1)=-Pointslr(:,1);end
        if(m(4)>m(3)), Pointsr(:,2)=-Pointsr(:,2); Pointslr(:,2)=-Pointslr(:,2);end
        [t,i]=sort(Evalues);
        Pointst=Pointsr(:,i);
        Pointslt=Pointslr(:,i);
    end
else
    if(options.rotate>1)
        % Eigen Vector align
        [Evalues, Evectors]=PCA(Pointst');
        
        % Do the Point transformation
        M=eye(3)/Evectors;
        Ml=zeros(4,4); Ml(1:3,1:3)=M; Ml(4,4)=1;
        Pointsr=transformvertices(Pointst,Ml);
		Pointslr=transformvertices(Pointslt,Ml);
		        
        % Heaviest axis flip
        R=sqrt(sum(Pointsr.^2,2));
        h=Pointsr>0;
        m(1)=sum(R(h(:,1))); m(2)=sum(R(~h(:,1)));
        m(3)=sum(R(h(:,2))); m(4)=sum(R(~h(:,2)));
        m(5)=sum(R(h(:,3))); m(6)=sum(R(~h(:,3)));
        if(m(2)>m(1)), Pointsr(:,1)=-Pointsr(:,1);  Pointslr(:,1)=-Pointslr(:,1);end
        if(m(4)>m(3)), Pointsr(:,2)=-Pointsr(:,2);  Pointslr(:,2)=-Pointslr(:,2);end
        if(m(6)>m(5)), Pointsr(:,3)=-Pointsr(:,3);  Pointslr(:,3)=-Pointslr(:,3);end
        [t,i]=sort(Evalues);
        Pointst= Pointsr(:,i);
        Pointslt= Pointslr(:,i);
    end
end
R1=mean(sqrt(sum(Pointst.^2,2)));
Pointsn=Pointst/R1;
Pointsln=Pointslt/R1;
    

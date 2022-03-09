function C=Features2CostMatrix_std(F1,F2,Points1,Points2,maxdist)
F1=F1'; %559*3375
F2=F2';%558*3375
C=zeros([size(F1,1) 1]);%559*558
for i=1:size(F1,1)
    P=F1(i,:)+eps+F2(i,:);%bsxfun(@plus,F1(:,i)+eps,F2(:,i)');
    M=F1(i,:)-F2(i,:);
%     M=bsxfun(@minus,F1(:,i),F2(:,i)');
    C_sub=((M.^2)./P);
    C(i)=sum(C_sub);
end
% P1=(Points1(:,1)-Points2(:,1)').^2;
% P2=(Points1(:,2)-Points2(:,2)').^2;
% P=P1+P2;
% C(P>maxdist^2)=inf;
function C=Features2CostMatrix(F1,F2,Points1,Points2,maxdist)
F1=F1'; %559*3375
F2=F2';%558*3375
C=zeros([size(F1,1) size(F2,1)]);%559*558
for i=1:size(F1,2)
    P=F1(:,i)+eps+F2(:,i)';%bsxfun(@plus,F1(:,i)+eps,F2(:,i)');
    M=F1(:,i)-F2(:,i)';
%     M=bsxfun(@minus,F1(:,i),F2(:,i)');
    C=C+((M.^2)./P);
end

% C=zeros([size(F1,1) size(F2,1) 3375]);%559*558
% for i=1:size(F1,1)
%     P=F1(i,:)+eps+F2;
%     M=F1(i,:)-F2;
%     value=((M.^2)./P);
% %     tic
% %     sum_value=sum(((M.^2)./P),2);
%     C(i,:,:)=value;
% end
% C=sum(C,3);

% P1=(Points1(:,1)-Points2(:,1)').^2;
% P2=(Points1(:,2)-Points2(:,2)').^2;
% P=P1+P2;
% C(P>maxdist^2)=inf;
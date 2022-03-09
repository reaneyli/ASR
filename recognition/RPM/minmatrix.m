function C=minmatrix(A)
[t,C]=min(A,[],2); C=C(:);
% C=zeros(size(A,1),1);
% for k=1:size(A,1);
%     [t,ind]=min(A(:));
%     [i,j]=ind2sub(size(A),ind);
%     A(i,:)=inf;
%     A(:,j)=inf;
%     C(i)=j;
% end
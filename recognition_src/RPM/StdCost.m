function [Cost]=StdCost(X_sub,X_tar,std)
 ntarcell=length(X_tar);
 nsubcell=length(X_sub);
std=std.^2; 
% calculate the assignment energy of each atlas point to all topre points
shape_prob=zeros(nsubcell,ntarcell);%row:assignment energy of one atlas point to all topre points
for i=1:3
    xyzdiff_atlas2topre=X_sub(:,i)-X_tar(:,i)';
    a=repmat(std(:,i),1,nsubcell);
    exponential_term=(xyzdiff_atlas2topre.^2-(a'));
    shape_prob=shape_prob+exponential_term;
end

Cost=shape_prob;
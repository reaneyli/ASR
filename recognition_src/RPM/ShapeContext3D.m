
function [Cost]=ShapeContext3D(Points1in,Points2in,Points1lin,Points2lin,options,Normals1,Normals2)
%
% [Matches,Cost]=ShapeContext(P1,P2,PL1,PL2,Options)
%
%  inputs,
%    P1 : Point List Nx2 or Nx3 describing the object contour in
%            the first dataset
%    P1 : Point List Mx2 or Mx3 describing the object contour in
%            the second dataset
%  (optional)
%    PL1 : Set equal to P1, or to a point list with more samples describing
%         the shame shape as P1
%    PL2 : Set equal to P2, or to a point list with more samples describing
%         the shame shape as P2
%    Options : A struct with options:
%      Options.r_bins : Number of Log Distance bins in the
%                           matching histogram(default 15)
%      Options.a_bins : Number of Angle bins in the
%                           matching histogram (default 15)
%      Options.rotate : Rotation Normalization of the points using
%                       eigen vector analyis 0,1 or 2.
%                       0: No correction (default), 
%                       1: Minimal align rotation
%                       2: Rotation with Heavyside flip.
%      Options.method : The matching method used:
%                       0 : Munkres One to one matching
%                       1 : Minimum One to multiple matching (default)
%      Options.maxdist : The maximum matching distance between normalized
%                        points (default 5). The Point-sets are normalized
%                        to have an average radius of one to their mean value.
%
%  outputs,
%    Matches : A 2xK list with the indices of the matched points
%    Cost : The histogram distance-cost for every match
%
%
% Example, 2D
%  % Add needed function paths
%   fnname='ShapeContext.m';
%   fndir=which(fnname); fndir=fndir(1:end-length(fnname));
%   addpath([fndir '/functions'])
%
%  % Read two images with bats inside
%   I1=imread('images\bat1.png');
%   I2=imread('images\bat2.png');
%
%  % Show the images
%   figure,
%   subplot(1,2,1), imshow(I1);
%   subplot(1,2,2), imshow(I2);
%
%  % Get the contours of the bats
%   [Lines1,Points1]=isocontour(I1,0.5);
%   [Lines2,Points2]=isocontour(I2,0.5);
%
%  % Do Shape Context matching, to find corresponding points
%   [Matches,Cost]=ShapeContext(Points1,Points2);
%
%  % Show the result
%   figure, hold on;
%   plot(Points1(:,2),Points1(:,1),'b.');
%   plot(Points2(:,2),Points2(:,1),'r.');
%   Cost=Cost./max(Cost(:));
%   P1=Points1(Matches(:,1),:);
%   P2=Points2(Matches(:,2),:);
%   for i=1:size(Matches,1)
%     C = Cost(i)*5;
%     plot([P1(i,2) P2(i,2)],[P1(i,1) P2(i,1)],'g','LineWidth',C);
%   end
%
% Example, 3D
%  % Load the 3D point datasets
%   load('images\sphere.mat');
%
%  % Do Shape Context matching, to find corresponding points
%   Matches=ShapeContext(Points1,Points2);
%
%  % Show the result
%   figure, hold on;
%   plot3(Points1(:,2),Points1(:,1),Points1(:,3),'b.');
%   plot3(Points2(:,2),Points2(:,1),Points2(:,3),'r.');
%   for i=1:size(Matches,1);
%     P1=Points1(Matches(i,1),:);
%     P2=Points2(Matches(i,2),:);
%     plot3([P1(2) P2(2)],[P1(1) P2(1)],[P1(3) P2(3)],'g');
%   end
%
%
% Example, 3D
%  % Load the 3D point datasets
%   load('testdata.mat');
%   Points1=FV1.vertices;
%   Points2=FV2.vertices;
%
%  % Do Shape Context matching, to find corresponding points
%   Matches=ShapeContext(Points1,Points2,FV1.faces,FV2.faces);
%
%  % Show the result
%   figure, hold on;
%   plot3(Points1(:,2),Points1(:,1),Points1(:,3),'b.');
%   plot3(Points2(:,2),Points2(:,1),Points2(:,3),'r.');
%   for i=1:size(Matches,1);
%     P1=Points1(Matches(i,1),:);
%     P2=Points2(Matches(i,2),:);
%     plot3([P1(2) P2(2)],[P1(1) P2(1)],[P1(3) P2(3)],'g');
%   end
%
% Function is written by D.Kroon University of Twente (February 2011)
% Process inputs
% defaultoptions=struct('r_max',2,'r_min',1e-9,'r_bins',15,'a_bins',15,'rotate',0,'method',1,'maxdist',5);
defaultoptions=struct('r_max',2,'r_min',1e-9,'r_bins',15,'a_bins',15,'rotate',0,'method',1,'maxdist',5);
options=defaultoptions;
% if(~exist('options','var')), options=defaultoptions;
% else
%     tags = fieldnames(defaultoptions);
%     for i=1:length(tags),
%         if(~isfield(options,tags{i})), options.(tags{i})=defaultoptions.(tags{i}); end
%     end
%     if(length(tags)~=length(fieldnames(options))),
%         warning('image_registration:unknownoption','unknown options found');
%     end
% end
% Check the Inputs
if(size(Points1in,2)<2||size(Points1in,2)>3), error('ShapeContext:inputs','Point list must be N x 2 or N x 3'); end
if(size(Points2in,2)<2||size(Points2in,2)>3), error('ShapeContext:inputs','Point list must be N x 2 or N x 3'); end
if(nargin<3), Points1lin=Points1in; Points2lin=Points2in; end
if(nargin<6), Normals1=[]; Normals2=[]; end
% Add needed Sub-folders to Path
fnname='ShapeContext.m';
fndir=which(fnname); fndir=fndir(1:end-length(fnname));
% addpath([fndir '/functions'])
% Normalize Translation, Scale and Orientation of Both Point clouds
Points1=Points1in;
Points1l=Points1in;
Points2=Points2in;
Points2l=Points2in;

% Small Alignment steps rotation 
if(options.rotate>0)
    [Points2,Points2l] = Horn_Align_ICP(Points1,Points2,Points2l);
end
% Determine the Distance-Polar histograms for every point.
% The histogram of a point describes the orientation and distance of
% all points relative to the point.
%
F1=getHistogramFeatures(Points1,Points1l,Normals1,options);
F2=getHistogramFeatures(Points2,Points2l,Normals2,options);

% Assign closest points. Process all points in once, or split the
% problem in point sub-volumes.
[Cost]=ShapeContextAligned(Points1,Points2,F1,F2,options);

function [Cost]=ShapeContextAligned(Points1,Points2,F1,F2,options)
% Calculate matrix on which every value give the matching cost between
% a certain point in dataset-1 and a certain point in dataset-2.
C=Features2CostMatrix(F1,F2,Points1,Points2,options.maxdist);
C(C==1e100)=inf;
Cost=C;

function C=Features2CostMatrix_Euclidean(F1,F2,Points1,Points2,maxdist)
F1=F1'; F2=F2';
C=zeros([size(F1,1) size(F2,1)]);
for i=1:size(F1,2)
    M=bsxfun(@minus,F1(:,i),F2(:,i)');
    C=C+M.^2;
end
% P1=bsxfun(@minus,Points1(:,1),Points2(:,1)').^2;
% P2=bsxfun(@minus,Points1(:,2),Points2(:,2)').^2;
% P=P1+P2;
% C(P>maxdist^2)=inf;

% function C=Features2CostMatrix(F1,F2,Points1,Points2,maxdist)
% F1=F1'; %559*3375
% F2=F2';%558*3375
% C=zeros([size(F1,1) size(F2,1)]);%559*558
% % for i=1:2:(size(F1,2)-1)
% %     P=F1(:,i)+eps+F2(:,i)';%bsxfun(@plus,F1(:,i)+eps,F2(:,i)');
% %     M=F1(:,i)-F2(:,i)';
% % %     M=bsxfun(@minus,F1(:,i),F2(:,i)');
% %     C=C+((M.^2)./P);
% %      P=F1(:,i+1)+eps+F2(:,i+1)';%bsxfun(@plus,F1(:,i)+eps,F2(:,i)');
% %     M=F1(:,i+1)-F2(:,i+1)';
% % %     M=bsxfun(@minus,F1(:,i),F2(:,i)');
% %     C=C+((M.^2)./P);
% % end
% for i=1:(size(F1,2))
%     P=F1(:,i)+eps+F2(:,i)';
%     M=F1(:,i)-F2(:,i)';
%     C=C+((M.^2)./P);
% end
% 
% % C=zeros([size(F1,1) size(F2,1) 3375]);%559*558
% % for i=1:size(F1,1)
% %     P=F1(i,:)+eps+F2;
% %     M=F1(i,:)-F2;
% %     value=((M.^2)./P);
% % %     tic
% % %     sum_value=sum(((M.^2)./P),2);
% %     C(i,:,:)=value;
% % end
% % C=sum(C,3);
% 
% P1=(Points1(:,1)-Points2(:,1)').^2;
% P2=(Points1(:,2)-Points2(:,2)').^2;
% P=P1+P2;
% C(P>maxdist^2)=inf;

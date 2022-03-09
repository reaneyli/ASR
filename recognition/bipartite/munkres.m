function [assignment,cost,dMat] = munkres(costMat)
% MUNKRES   Munkres Assign Algorithm 
%
% [ASSIGN,COST] = munkres(COSTMAT) returns the optimal assignment in ASSIGN
% with the minimum COST based on the assignment problem represented by the
% COSTMAT, where the (i,j)th element represents the cost to assign the jth
% job to the ith worker.
%

% This is vectorized implementation of the algorithm. It is the fastest
% among all Matlab implementations of the algorithm.

% Examples
% Example 1: a 5 x 5 example
%{
[assignment,cost] = munkres(magic(5));
[assignedrows,dum]=find(assignment);
disp(assignedrows'); % 3 2 1 5 4
disp(cost); %15
%}
% Example 2: 400 x 400 random data
%{
n=400;
A=rand(n);
tic
[a,b]=munkres(A);
toc                 % about 6 seconds 
%}

% Reference:
% "Munkres' Assignment Algorithm, Modified for Rectangular Matrices", 
% http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html

% version 1.0 by Yi Cao at Cranfield University on 17th June 2008

assignment = false(size(costMat));
cost = 0;

costMat(costMat~=costMat)=Inf;
validMat = costMat<Inf;
validCol = any(validMat);
validRow = any(validMat,2);

nRows = sum(validRow);
nCols = sum(validCol);
n = max(nRows,nCols);
if ~n
    return
end
    
dMat = zeros(n);
dMat(1:nRows,1:nCols) = costMat(validRow,validCol);  %�������

%*************************************************
% Munkres' Assignment Algorithm starts here
%*************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   STEP 1: Subtract the row minimum from each row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 dMat = bsxfun(@minus, dMat, min(dMat,[],2));   %��ȥÿ�е���Сֵ


%**************************************************************************  
%   STEP 2: Find a zero of dMat. If there are no starred zeros in its
%           column or row start the zero. Repeat for each zero
%**************************************************************************
zP = ~dMat;  %ÿ����Сֵ��λ�ø�ֵΪ1�������ط���Ϊ0
starZ = false(n); %�ؽ�n*n ��0����
%zP�����а�������ȡ��һ�γ��ֵ���ֵλ��Ϊ1����˼����������ֶ�����ͬһ�ж�Ϊ��Сֵ����ôȡ��һ�γ��ֵ���һ����Ϊ��Сֵ,����starZ��ֵΪ1��
%ʣ��δ����ֵ���н�ȫ��0
%����Ϊ��Сֵ���ҵ�һ�γ��ֵ�λ�ö���������starZ��
while any(zP(:))
    [r,c]=find(zP,1);%�ҳ���һ������Ԫ�ص�λ��
    starZ(r,c)=true;%starZ�и�λ��Ԫ�ظ�Ϊ1
    zP(r,:)=false;%��λ��Ԫ�ص��к��ж�Ϊ0
    zP(:,c)=false;
end

while 1
%**************************************************************************
%   STEP 3: Cover each column with a starred zero. If all the columns are
%           covered then the matching is maximum
%**************************************************************************
    primeZ = false(n);
    coverColumn = any(starZ);%�ж�ÿһ���Ƿ��з���Ԫ�أ���Щ������������Сֵ��λ������ǰ��Щ���ظ�
    if ~any(~coverColumn)
        break
    end
    coverRow = false(n,1);%����n*1�ľ���
    while 1
        %**************************************************************************
        %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
        %           zero in the row containing this primed zero, Go to Step 5.  
        %           Otherwise, cover this row and uncover the column containing 
        %           the starred zero. Continue in this manner until there are no 
        %           uncovered zeros left. Save the smallest uncovered value and 
        %           Go to Step 6.
        %**************************************************************************
        zP(:) = false;
        zP(~coverRow,~coverColumn) = ~dMat(~coverRow,~coverColumn);%�ҳ�ZP�ж��г�����Сֵ������δ��ѡ���λ�ã����ڸ�λ����Ȼ����ֵΪ1��
        Step = 6;
		%ʣ���е��ظ�����Сֵ������ȡ�ʼ���ֵ�����Ϊ��Сֵ
        while any(any(zP(~coverRow,~coverColumn)))
            [uZr,uZc] = find(zP,1);%�ҳ�����λ��
            primeZ(uZr,uZc) = true;%����Щλ�ø���primeZ��
            stz = starZ(uZr,:);%�������ԭ���ҵ�����Сֵ���е�����Ԫ��
            if ~any(stz)
                Step = 5;
                break;
            end
            coverRow(uZr) = true;%�洢�������ظ���Ԫ�ص���λ����Ϣ��0->1
            coverColumn(stz) = false;%���ҵ��ĸ���Ԫ��תΪ0����˼�����ҳ��ʼ���ظ��е��У��洢�������ظ�Ԫ�ص���
            zP(uZr,:) = false;%ȥ���ҵ����ظ�Ԫ�����ڵ���
            zP(~coverRow,stz) = ~dMat(~coverRow,stz); %��dMat��stz��ʣ���Ԫ�ظ�ֵ��ZP�У�����Ƿ�stz���Ƿ���0Ԫ�أ����ƹ�ȥ���1
        end
        if Step == 6
            % *************************************************************************
            % STEP 6: Add the minimum uncovered value to every element of each covered
            %         row, and subtract it from every element of each uncovered column.
            %         Return to Step 4 without altering any stars, primes, or covered lines.
            %**************************************************************************
            M=dMat(~coverRow,~coverColumn);
            m1=dMat(coverRow,coverColumn);
            minval=min(min(M));
            if minval==inf
                return
            end
            dMat(coverRow,coverColumn)=dMat(coverRow,coverColumn)+minval;
            dMat(~coverRow,~coverColumn)=M-minval;
            dMat=dMat;
        else
            break
        end
    end
    %**************************************************************************
    % STEP 5:
    %  Construct a series of alternating primed and starred zeros as
    %  follows:
    %  Let Z0 represent the uncovered primed zero found in Step 4.
    %  Let Z1 denote the starred zero in the column of Z0 (if any).
    %  Let Z2 denote the primed zero in the row of Z1 (there will always
    %  be one).  Continue until the series terminates at a primed zero
    %  that has no starred zero in its column.  Unstar each starred
    %  zero of the series, star each primed zero of the series, erase
    %  all primes and uncover every line in the matrix.  Return to Step 3.
    %**************************************************************************
    rowZ1 = starZ(:,uZc);
    starZ(uZr,uZc)=true;
    while any(rowZ1)
        starZ(rowZ1,uZc)=false;
        uZc = primeZ(rowZ1,:);
        uZr = rowZ1;
        rowZ1 = starZ(:,uZc);
        starZ(uZr,uZc)=true;
    end
end

% Cost of assignment
assignment(validRow,validCol) = starZ(1:nRows,1:nCols);
cost = sum(costMat(assignment));
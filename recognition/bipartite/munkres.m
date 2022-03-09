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
dMat(1:nRows,1:nCols) = costMat(validRow,validCol);  %补零操作

%*************************************************
% Munkres' Assignment Algorithm starts here
%*************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   STEP 1: Subtract the row minimum from each row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 dMat = bsxfun(@minus, dMat, min(dMat,[],2));   %减去每行的最小值


%**************************************************************************  
%   STEP 2: Find a zero of dMat. If there are no starred zeros in its
%           column or row start the zero. Repeat for each zero
%**************************************************************************
zP = ~dMat;  %每行最小值的位置赋值为1，其他地方都为0
starZ = false(n); %重建n*n 的0矩阵
%zP矩阵中按列依次取第一次出现的最值位置为1，意思就是如果出现多行中同一列都为最小值，那么取第一次出现的那一列作为最小值,并在starZ赋值为1，
%剩下未被赋值的列将全是0
%所有为最小值，且第一次出现的位置都被保存在starZ中
while any(zP(:))
    [r,c]=find(zP,1);%找出第一个非零元素的位置
    starZ(r,c)=true;%starZ中该位置元素赋为1
    zP(r,:)=false;%该位置元素的行和列都为0
    zP(:,c)=false;
end

while 1
%**************************************************************************
%   STEP 3: Cover each column with a starred zero. If all the columns are
%           covered then the matching is maximum
%**************************************************************************
    primeZ = false(n);
    coverColumn = any(starZ);%判断每一列是否有非零元素，这些列则是由于最小值的位置与先前有些行重复
    if ~any(~coverColumn)
        break
    end
    coverRow = false(n,1);%定义n*1的矩阵
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
        zP(~coverRow,~coverColumn) = ~dMat(~coverRow,~coverColumn);%找出ZP中多行出现最小值，但是未被选择的位置，现在该位置依然被赋值为1，
        Step = 6;
		%剩下中的重复的最小值的行中取最开始出现的列作为最小值
        while any(any(zP(~coverRow,~coverColumn)))
            [uZr,uZc] = find(zP,1);%找出坐标位置
            primeZ(uZr,uZc) = true;%把这些位置赋到primeZ中
            stz = starZ(uZr,:);%并且提出原来找到的最小值该行的所有元素
            if ~any(stz)
                Step = 5;
                break;
            end
            coverRow(uZr) = true;%存储了所有重复的元素的行位置信息，0->1
            coverColumn(stz) = false;%将找到的该列元素转为0，意思则是找出最开始有重复行的列，存储了所有重复元素的列
            zP(uZr,:) = false;%去除找到的重复元素所在的行
            zP(~coverRow,stz) = ~dMat(~coverRow,stz); %把dMat中stz中剩余的元素赋值到ZP中，检测是否stz中是否还有0元素，复制过去变成1
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

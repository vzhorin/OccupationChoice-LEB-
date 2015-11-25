% Function that merges two distributions and its supports into one of 
% equal dimension
%1:no credit : created by resdistnc
%2:credit: created by resdistc


function [sp,pr] = mdist(sp1,sp2,pr1,pr2,F,pts,bds)

% THIS IS WEIRD: WE USE INPUTS NOT DEFINED AND ALSO TOO MANY INPUT ARGUMENTS !!

%nargin=5  actually

if (nargin == 5)
    pts     = length(sp1);% 100 points
    
end;
sp      = union(sp1,sp2);% size (1,200)

if (nargin ~= 7) % we are not in this case
     smax    = length(sp);
     sp      = linspace(sp(1),sp(smax),pts);
else % we are in this case
    
     sp      = linspace(bds(1),bds(2),pts);
 
 end;
n1      = hist(sp1,sp);% noncredit
n2      = hist(sp2,sp);% credit
nc1     = cumsum(n1);% noncredit
nc2     = cumsum(n2);% credit
pr      = zeros(1,pts);% 200 points to store the pdf
jp1=1;
jp2=1;
for i = 1:pts %100 points
    jn1=nc1(i);
    jn2=nc2(i);
    pr(i) = sum(pr1(jp1:jn1))*(1-F);% noncredit
    pr(i) = pr(i) + sum(pr2(jp2:jn2))*F;% credit
    jp1=jn1+1;
    jp2=jn2+1;
end;

%size(sp) (1,100)
%size(pr) (1,100)






function [x,z_max]=simplex(z,A,b)

%The objective function's coefficients are a row vector z; the left-hand
%coeffecients of the main constraints are A; the right-hand of the main
%constrants are a column vector b'. We assume main variables x1,..,xn are
%nonnegative. We also assume that the constraint matrix A is m-by-n.

%The output is a solution vector x in R^n and the resulting objective
%function.

%Maybe we should introduce a tolerance for this function to improve run
%time.

%Initialize augmented matrix to form initial simplex tableau
sz=size(A);
m=sz(1,1);
n=sz(1,2);
I=eye(m,n);
Z=zeros(m,1);
b=b';
%Initialize initial simplex tableau (ITS) in augmented matrix form
M0=[Z A I];
M1=cat(2,M0,b);
M2=[-1 z zeros(1,m) 0];
M=cat(1,M1,M2);

%Identify pivot of ITS
obj=M(end,:);
[value,col]=max(obj);
p=M(:,col);
p(m+1)=[];
ratio=p./b;
r0=min(ratio);
[row]=find(ratio == r0);
pivot=M(row,col);

%Conduct first simplex pivot
M(row,:)=(1/pivot)*M(row,:);
for i=row+1:m+1
    M(i,:)=M(i,:)-M(i,col).*M(row,:);
end

while 1
    %Examine objective function
    obj=M(end,:);
    %Identify next pivot column, y and the corresponding coefficient rho
    y=find(obj==max(obj));
    rho=max(obj);
    %Compute vector of ratios
    chi=zeros(1,m);
    for i=1:m
        chi(i)=(M(i,end))/rho;
    end
    %Find pivot row, w
    w=find(chi==min(chi));
    pivot=M(w,y);
    
    %Pivot to next simplex tableau
    M(w,:)=(1/pivot)*M(w,:);
    for i=w+1:m+1
        M(i,:)=M(i,:)-M(i,y).*M(w,:);
    end
    
    %Stop iteration when objective function has all nonpositive
    %coefficients
    obj=M(end,:);
    if max(obj)==0
        break
    end
end
%Once optimal tableau is reached, use objective function to write solution
z_max=-M(end,end);
obj=M(end,:);
x=zeros(1,2*n);
sigma1=find(obj==0);
sigma=sigma1-1;
for i=1:length(sigma)
    x(sigma(i))=M(i,end);
end
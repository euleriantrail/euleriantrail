function [x]=conjugate_gradient(A,b)

%From Trefethen

%This functionsolves for unknown vectors in linear equations of the form 
%Ax=b, where A is n-by-n, real, symmetric, psd, and b is a vector in R^n.

%Initialize
sz=size(A);
n=sz(1,1);
x=zeros(n,1);
r=b;
p=r;

for i=1:100
    %Step length
    alpha=(r'*r)/(p'*A*p);
    %Approximate solution
    x=x+alpha*p;
    %Residual
    r=r-alpha*A*p;
    %Improvement this step
    beta=(r'*r)/(r'*r);
    %Search direction
    p=r+beta*p;
end
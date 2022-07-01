function K=kerf(x1,x2,mu)
% Compute: K=K(A,B);
% Gaussian Kernel: k(x,y)=exp^(-mu||x-y||^2)
[m1,n1]=size(x1);
[m2,n2]=size(x2);
if n1~=n2
    disp('kernel product error: n1 != n2');
    return;
end
% if m2>800
%     x2=x2(1:100,:);
%     m2=size(x2,1);
% end
x11=sum(x1.^2,2)*ones(1,m2);
x22=sum(x2.^2,2)*ones(1,m1);
K=x11+x22'-2*(x1*x2');
K=exp(-mu*K);
end
%% This function is used for HGMM, but need to selected
function [a,b,c]= DualSolver(Y,K,dK,lambda)
opt=optimoptions('quadprog','Display','off');
m=length(Y);
Yt=Y';
yky=Y.*K.*Yt;  %YX'XY
ydy=Y.*dK.*Yt; %YDeltaY/lambda
H=yky+Y.*(dK*dK').*Yt/lambda;
e=ones(m,1);
warning off all;
a = quadprog(H,-e,[],[],Yt,0,zeros(m,1),[],[],opt);
c=Y.*dK'.*Yt*a/lambda;
tol=1e-6;
ind=find(a>tol);
b=Y(ind,:)-Y(ind,:).*(yky(ind,:)*a+ydy(ind,:)*c);
b=mean(b);
end

%% This function is used for SGMM, but need to selected
function [a,b,c]= DualSolver(Y,K,dK,lambda,C)
opt=optimoptions('quadprog','Display','off');
m=length(Y);
Yt=Y';
yky=Y.*K.*Yt;  %YX'XY
ydy=Y.*dK.*Yt; %YDeltaY/lambda
H=yky+Y.*(dK*dK').*Yt/lambda;
e=ones(m,1);
warning off all;
a = quadprog(H,-e,[],[],Yt,0,zeros(m,1),C*e,[],opt);
c=Y.*dK'.*Yt*a/lambda;
tol=1e-6;
ind=find(a>tol & a<C-tol);
if ~isempty(ind)
    b=Y(ind,:)-Y(ind,:).*(yky(ind,:)*a+ydy(ind,:)*c);
    b=mean(b);
else
    wX=1-yky*a-ydy*c*lambda;  
    bl1=max(wX(a<=tol & Y==1));
    bu1=min(-wX(a<=tol & Y~=1));
    bu2=min(wX(a>=c-tol & Y==1));
    bl2=max(-wX(a>=c-tol & Y==-1));
    b=(max([bl1,bl2])+min([bu1,bu2]))/2;
    if isempty(b)
        b=0;
    end
end
end


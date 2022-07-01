function pY = HGMM(testX,X,Y,lambda,mu,dmu)
%MSVMMAIN 
% mu for nonlinear case, and =0 is linear case
% dmu for delta function, delta(x_i,x_j)=exp{-dmu*||x_i-x_j||^2}
if mu~=0
    K=kerf(X,X,mu);
else
    K=X*X';
end
dK=kerf(X,X,dmu);


[a,b,c]= DualSolver(Y,K,dK,lambda);
pY= Predict(testX,X,Y,a,b,c,mu,dmu);
end


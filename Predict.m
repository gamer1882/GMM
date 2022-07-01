function pY = Predict(testX,X,Y,a,b,c,mu,dmu)
%SVMPREDICT 
pY=zeros(size(testX,1),1);
tol=1e-6;
ind1=find(a>tol);
ind2=find(c>tol);
if mu~=0
    pK=kerf(testX,X(ind1,:),mu);
else
    pK=testX*X(ind1,:)';
end
pdK=kerf(testX,X(ind2,:),dmu);


fX=pK*(a(ind1).*Y(ind1))+b+pdK*(c(ind2).*Y(ind2));
pY(fX>eps)=1;
pY(fX<-eps)=-1;
end


function x=noisyIRLS(y,D,p,tr,lam,maxiter)

K=size(D,2);
W=eye(K);
x=zeros(K,1);
DTD=D'*D;
DTy=D'*y;
i=1;
while 1
   x0=(lam*W+DTD)\DTy;
   W=diag(abs(x0).^(2-p)+eps);
   if norm(x0-x) <= tr || i > maxiter
       break
   end
   x=x0;
   i=i+1;
end
x=x0;
return
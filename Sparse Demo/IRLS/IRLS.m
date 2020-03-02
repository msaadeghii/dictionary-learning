function x=IRLS(y,D,p,tr,maxiter)

K=size(D,2);
W=eye(K);
x=zeros(K,1);
i=1;
while 1
   P1=D*W*D';
   P2=W*D';
   x0=P2*pinv(P1)*y;
   W=diag(abs(x0).^(2-p)+eps);
   if norm(x0-x) <= tr || i > maxiter
       break
   end
   x=x0;
   i=i+1;
end
x=x0;
return
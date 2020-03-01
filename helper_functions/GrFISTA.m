function X=GrFISTA(B,A,X,lam,maxit,map)
% Group (batch) fast iteative shrinkage thresholding algorithm (FISTA) for sparse representation
AtA=A'*A;
AtB=A'*B;
t = 1/(2*eigs(AtA,1));
tr=1e-5;
L=size(B,2);
Y=X;
thetp=1;
it=0;
er=1;

while er && (it<=maxit)
    Xp=X;
    Z=Y-t*AtA*Y+t*AtB;
    Zb=abs(Z);Si=Zb>=(lam*t);
    X=(Z-lam*t*sign(Z)).*Si;
    thet=(1+sqrt(1+4*thetp^2))/2;
    bet=(thetp-1)/thet;
    Xdiff=(X-Xp);
    Y=X+bet*Xdiff;
    thetp=thet;
    it=it+1;
    er= sum(sum((Xp-X).^2))>tr;
end

if map == 1,
% debiasing
  
    for l = 1:L,
        x = (abs(X(:,l))>1e-5);
        M = diag(x);
        ANew = A*M;
        AS = ANew(:,x);
        X(x,l) =(AS'*AS+.0001*eye(sum(x)))\(AS'*B(:,l));
    end    
end

return
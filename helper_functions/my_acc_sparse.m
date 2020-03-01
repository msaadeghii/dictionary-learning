function X=my_acc_sparse(Y,D,X,lam,tau,w)


diff=inf;
mu=1/svds(D,1)^2;
Xo=X;

while diff > tau
    
    T=X+w*(X-Xo);
    Xo=X;
    X=my_softh(T-mu*D'*(D*T-Y), lam*mu); 
    diff=norm(X-Xo,'fro')/norm(Xo,'fro');
    
end


end
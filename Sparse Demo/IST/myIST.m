function x=myIST(y,D,lam,er)
    m=size(D,2);
    x=pinv(D)*y;
    DTy=D'*y;
    DTD=D'*D;
    xo=ones(m,1);
    d=eigs(DTD);
    c=d(1)+0.05;

    while norm(x-xo)>er
        xo=x;
        x=softh(x+1/c*(DTy-DTD*x),lam/c);   
    end;
    ind=abs(x)>1e-5;
    x=zeros(m,1);
    x(ind)=pinv(D(:,ind))*y;
end
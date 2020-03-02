function [x] = conjgrad(A,b,x)
    r=b-A*x;
    p=r;
    rsold=r'*r;
 
    for i=1:size(A,1)
        Ap=A*p;
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        if sqrt(rsnew)<1e-10
              break;
        end
        p=r+rsnew/rsold*p;
        rsold=rsnew;
    end
end
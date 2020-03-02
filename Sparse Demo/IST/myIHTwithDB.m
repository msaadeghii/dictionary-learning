function x=myIHTwithDB(y,D,s,lam,er)
    m=size(D,2);
    x=pinv(D)*y;
    DTy=D'*y;
    DTD=D'*D;
    xo=ones(m,1);

    while norm(x-xo)>er
        xo=x;
        [~,ind]=selectsl(x+lam*(DTy-DTD*x),s);
        x=zeros(m,1);
        x(ind)=pinv(D(:,ind))*y;
    end
    
    function [yo ind]=selectsl(y,s)
        yab=abs(y);
        [a b]=sort(yab,'descend');
        yo=y;
        yo(yab<a(s))=0;
        ind=b(1:s);
       
    end
end
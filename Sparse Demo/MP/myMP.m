function x=myMP(y,D,s)
    m=size(D,2);
    x=zeros(m,1);
    r=y;
    for i=1:s
        c=abs(D'*r);
        ind=find(c==max(c),1);
        x(ind)=x(ind)+D(:,ind)'*r;
        r=r-D(:,ind)*D(:,ind)'*r;
    end;
end
function x=myOMP(y,D,s)
    m=size(D,2);
    r=y;
    S=[];
    for i=1:s
        c=abs(D'*r);
        ind=find(c==max(c));
        S=sort([S,ind(1)]);
        xi=pinv(D(:,S))*y;
        r=y-D(:,S)*xi;    
    end;
    x=zeros(m,1);
    x(S)=xi;
end
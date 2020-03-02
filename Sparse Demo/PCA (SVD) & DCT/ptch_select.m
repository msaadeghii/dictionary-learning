function ptch=ptch_select(im,K)
[M N]=size(im);
p=0;bb=8;stp=1;
for i=1:stp:M-bb
    for j=1:stp:N-bb
        p=p+1;
        y=im(i:i+bb-1,j:j+bb-1);  
        Sn(p)=std(y(:));
    end
end
bins=200;
h2=hist(Sn,bins);
t2=min(Sn):(max(Sn)-min(Sn))/(bins-1):max(Sn);
m=K*mean(h2)/5;
T=m./(1+h2);
ind=T>1;
T(ind)=1;
p=0;    
for i=1:stp:M-bb
    for j=1:stp:N-bb
        y=reshape(im(i:i+bb-1,j:j+bb-1),bb^2,1);  
        Df=t2-std(y(:));
        [~,ind]=min(abs(Df));
        rnd=rand(1);
        if rnd>1-T(ind)^2.6
        p=p+1;
        ptch(:,p)=y;
        end
    end
end
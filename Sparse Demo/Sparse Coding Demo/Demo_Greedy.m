% Figures 3.5 and 3.6 of Elad's book. This code is from M. Elad's homepage.
% =========================================
% In this program we demonstrate the MP, OMP, WMP, and Thresh
% algorithm, by running them on a set of test signals and checking 
% wether they provide the desired outcome
clc
n=20; m=50; Smax=10; Exper=100; 

A=randn(n,m);
W=sqrt(diag(A'*A));
for k=1:1:m, 
    A(:,k)=A(:,k)/W(k); 
end; 

Er2=zeros(Smax,Exper,4); 
ErS=zeros(Smax,Exper,4); 

for S=1:1:Smax,
    for experiment=1:1:Exper
       
        % Generate a test signal of cardinality S
        x=zeros(m,1);
        pos=randperm(m);
        x(pos(1:S))=sign(randn(S,1)).*(1+rand(S,1));
        b=A*x; 

        % Apply LS-MP
        tic;
        thrLSMP=1e-4;
        r=b;
        SS=[];
        while r'*r>thrLSMP,
            Z=zeros(m,1);
            for jj=1:1:m
                SStemp=[SS,jj]; 
                rtemp=b-A(:,SStemp)*pinv(A(:,SStemp))*b;
                Z(jj)=rtemp'*rtemp;
            end;
            posZ=find(Z==min(Z),1);
            SS=sort([SS,posZ(1)]);
            r=b-A(:,SS)*pinv(A(:,SS))*b;    
        end;
        xLSMP=zeros(m,1);
        xLSMP(SS)=pinv(A(:,SS))*b;
        Er2(S,experiment,1)=mean((xLSMP-x).^2)/mean(x.^2);
        ErS(S,experiment,1)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));
        t_LSOMP=toc;
        % Apply OMP
        tic;
        thrOMP=1e-4;
        r=b;
        SS=[];
        while r'*r>thrOMP,
            Z=abs(A'*r);
            posZ=find(Z==max(Z));
            SS=sort([SS,posZ(1)]);
            r=b-A(:,SS)*pinv(A(:,SS))*b;    
        end;
        xOMP=zeros(m,1);
        xOMP(SS)=pinv(A(:,SS))*b;
        Er2(S,experiment,2)=mean((xOMP-x).^2)/mean(x.^2);
        ErS(S,experiment,2)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));
        t_OMP=toc;
        % Apply MP
        tic;
        thrMP=1e-4;
        r=b;
        xMP=zeros(m,1);
        while r'*r>thrMP,
            Z=abs(A'*r);
            posZ=find(Z==max(Z),1);
            xMP(posZ)=xMP(posZ)+A(:,posZ)'*r;
            r=r-A(:,posZ)*A(:,posZ)'*r;
        end;
        SS=find(abs(xMP)>1e-8)';
        Er2(S,experiment,3)=mean((xMP-x).^2)/mean(x.^2);
        ErS(S,experiment,3)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));
        t_MP=toc;
        % Apply WMP
        tic;
        thrWMP=1e-4; t=0.5; 
        r=b;
        xWMP=zeros(m,1);
        while r'*r>thrMP,
            Z=abs(A'*r);
            posZ=find(Z>=t*sqrt(r'*r),1);
            if isempty(posZ)
                posZ=find(Z==max(Z),1);
            end;
            xWMP(posZ)=xWMP(posZ)+A(:,posZ)'*r;
            r=r-A(:,posZ)*A(:,posZ)'*r;
        end;
        SS=find(abs(xWMP)>1e-8)';
        Er2(S,experiment,4)=mean((xWMP-x).^2)/mean(x.^2);
        ErS(S,experiment,4)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));
        t_WMP=toc;
        % Apply Thr 
        tic;
        thrTH=1e-4;
        Z=A'*b;
        [Za,posZ]=sort(abs(Z),'descend');
        % in case the support is known:
        % SS=posZ(1:S)'; xTH=zeros(m,1); 
        % xTH(SS)=pinv(A(:,SS))*b; 
        r=b;
        SS=[];
        xTH=zeros(m,1);
        while r'*r>thrTH,
            SS=[SS,posZ(length(SS)+1)];
            xTH=zeros(m,1);
            xTH(SS)=pinv(A(:,SS))*b;
            r=b-A(:,SS)*xTH(SS);
        end;
        Er2(S,experiment,5)=min(mean((xTH-x).^2)/mean(x.^2),1);
        ErS(S,experiment,5)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));
        t_Thr=toc;
        
    end;
    
    % displaying the results for the set of tests
    disp([mean(Er2(S,:,1)),mean(Er2(S,:,2)),...
            mean(Er2(S,:,3)),mean(Er2(S,:,4)),mean(Er2(S,:,5)) ... 
             mean(ErS(S,:,1)),mean(ErS(S,:,2)),...
             mean(ErS(S,:,3)),mean(ErS(S,:,4)),mean(ErS(S,:,5))]); 
    
end;

figure(1); clf; 
h=plot(1:1:Smax,mean(Er2(:,:,1),2),'b'); hold on; 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,2),2),'r'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,3),2),'g'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,4),2),'c'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,5),2),'m'); 
set(h,'LineWidth',2);
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Average and Relative L_2-Error'); set(h,'FontSize',14);
set(gca,'FontSize',14);
h=legend({'LS-OMP','OMP','MP','Weak-MP (t=0.5)','Thresholding'},2);
axis([0 Smax 0 1]);grid on
% print -depsc2 Chapter_03_GreedyDemoL2.eps

figure(2); clf; 
h=plot(1:1:Smax,mean(ErS(:,:,1),2),'b'); hold on; 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(ErS(:,:,2),2),'r'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(ErS(:,:,3),2),'g'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(ErS(:,:,4),2),'c'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(ErS(:,:,5),2),'m'); 
set(h,'LineWidth',2);
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Probability of Error in Support'); set(h,'FontSize',14);
set(gca,'FontSize',14);
h=legend({'LS-OMP','OMP','MP','Weak-MP (t=0.5)','Thresholding'},2);
grid on
title(['LS-OMP: ',num2str(t_LSOMP),'-- OMP: ',num2str(t_OMP),'-- MP: ',num2str(t_MP),'-- WMP: ',num2str(t_WMP),'-- Thr: ',num2str(t_Thr)]);
% print -depsc2 Chapter_03_GreedyDemoSupp.eps





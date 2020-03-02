% Figures 3.8 and 3.9 of Elad's book. This code is from M. Elad's homepage.
% =========================================
% In this program we demonstrate the OMP, IRLS and BP 
% algorithms, by running them on a set of test signals and checking 
% wether they provide the desired outcome
%=========================================

n=20; m=50; Smax=10; Exper=10; 

A=randn(n,m);
W=sqrt(diag(A'*A));
for k=1:1:m, 
    A(:,k)=A(:,k)/W(k); 
end; 

Er2=zeros(Smax,Exper,4); 
ErS=zeros(Smax,Exper,4); 
OMPtime=0; 
BPtime=0;
IRLStime=0;

for S=1:1:Smax,
    for experiment=1:1:Exper
       
        % Generate a test signal of cardinality S
        x=zeros(m,1);
        pos=randperm(m);
        x(pos(1:S))=sign(randn(S,1)).*(1+rand(S,1));
        b=A*x; 

        % Apply OMP
        thrOMP=1e-4;
        tic; 
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
        OMPtime=OMPtime+double(toc); 
        Er2(S,experiment,1)=mean((xOMP-x).^2)/mean(x.^2);
        ErS(S,experiment,1)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));
        
        % Apply IRLS
        tic;
        xIRLS=ones(m,1); % initialization
        XX=eye(m); 
        p=0.5; % refers to L1 (L0 is obtained by p=1)
        for k=1:1:100,
            XX2=XX*XX;
            xIRLS=XX2*A'*pinv(A*XX2*A')*b;
            XX=diag(abs(xIRLS).^p);
        end;
        SS=find(abs(xIRLS)>1e-4)'; 
        IRLStime=IRLStime+double(toc); 
        Er2(S,experiment,2)=mean((xIRLS-x).^2)/mean(x.^2);
        ErS(S,experiment,2)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));
        
        % BP using L1 by Matlab
        V=ones(2*m,1);
        tic;
        xBP=linprog(V,[],[],[A, -A],b,0*V,V*100);
        xBP=xBP(1:m)-xBP(m+1:end);
        SS=find(abs(xBP)>1e-4)'; 
        BPtime=BPtime+double(toc); 
        Er2(S,experiment,3)=mean((xBP-x).^2)/mean(x.^2);
        ErS(S,experiment,3)=(max(S,length(SS))-...
                              length(intersect(SS,pos(1:S))))/max(S,length(SS));

    end;
    
    % displaying the results for the set of tests
    disp([mean(Er2(S,:,1)),mean(Er2(S,:,2)),...
            mean(Er2(S,:,3)), ... 
             mean(ErS(S,:,1)),mean(ErS(S,:,2)),...
             mean(ErS(S,:,3))]); 
    
end;

figure(1); clf; 
h=plot(1:1:Smax,mean(Er2(:,:,1),2),'k'); hold on; 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,2),2),'k-.'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,3),2),'ko'); 
set(h,'MarkerFaceColor','k');
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Average and Relative L_2-Error'); set(h,'FontSize',14);
set(gca,'FontSize',14);
h=legend({'OMP','IRLS','BP by Linear Prog.'},2);
grid on
% print -depsc2 Chapter_03_RelaxationDemoL2.eps

figure(2); clf; 
h=plot(1:1:Smax,mean(ErS(:,:,1),2),'k'); hold on; 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(ErS(:,:,2),2),'k-.'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(ErS(:,:,3),2),'ko'); 
set(h,'MarkerFaceColor','k');
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Probability of Error in Support'); set(h,'FontSize',14);
set(gca,'FontSize',14);
h=legend({'OMP','IRLS','BP by Linear Prog.'},2);
grid on
% print -depsc2 Chapter_03_RelaxationDemoSupp.eps





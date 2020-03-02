% Figures - 5.6, 5.7, 5.8, and 5.9  of Elad's book. This code is from M. Elad's homepage.
% =========================================
% This program presents the IRLS solution to the Q problem. 
% We demonstrate several things here:
%       1. Effect of lambda and the optimal value
%       2. Kinds of solutions obtained for various lambda
%       3. The convergence of the IRLS
%       4. The solution path (LARS style)

% Creating the data to test on
n=100; m=200; S=4; sigma=0.1;

A=randn(n,m);
W=sqrt(diag(A'*A));
for k=1:1:m,
    A(:,k)=A(:,k)/W(k);
end;

x0=zeros(m,1);
pos=randperm(m);
x0(pos(1:S))=sign(randn(S,1)).*(1+rand(S,1));
b=A*x0+randn(n,1)*sigma;

save temp n m A b x0 S sigma pos

% The IRLS algorithm for varying lambda
Xsol=zeros(m,101);
Err=zeros(101,1);
Res=zeros(101,1);
lambda=0.0001*2.^(0:0.15:15);
AA=A'*A;
Ab=A'*b;
for ll=1:1:101
    xout=zeros(m,1);
    XX=eye(m);
    for k=1:1:15
        xout=inv(2*lambda(ll)*XX+AA)*Ab;
        XX=diag(1./(abs(xout)+1e-10));
        % disp(norm(xout-x0)/norm(x0));
    end;
    Xsol(:,ll)=xout; 
    Err(ll)=(xout-x0)'*(xout-x0)/(x0'*x0);
    Res(ll)=(A*xout-b)'*(A*xout-b)/sigma^2/n; 
    disp(Err(ll));
end;

figure(1); clf; 
h=semilogx(lambda,Err,'k');
set(h,'LineWidth',2); 
hold on;
h=semilogx(lambda,abs(log10(Res)),'k--');
set(h,'LineWidth',2); 
h=xlabel('\lambda'); 
set(h,'FontSize',14);
h=ylabel('Normalized L_2-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
grid on;
axis([min(lambda) max(lambda) 0 1]); 
% print -depsc2 Chapter_05_IRLS_lambdaCurve.eps

figure(2); clf;
posL=find(Err==min(Err),1); 
subplot(3,1,1); h=plot(Xsol(:,posL-22),'b'); 
set(h,'LineWidth',2);
hold on; h=plot(pos(1:S),x0(pos(1:S)),'ro'); 
set(h,'MarkerFaceColor','r');
set(gca,'FontSize',14);
h=text(50,0.65,['\lambda=',num2str(lambda(posL-22))]);
set(h,'FontSize',14);
subplot(3,1,2); h=plot(Xsol(:,posL)); 
set(h,'LineWidth',2);
hold on; h=plot(pos(1:S),x0(pos(1:S)),'ro'); 
set(h,'MarkerFaceColor','r');
set(gca,'FontSize',14);
h=text(50,0.65,['\lambda=',num2str(lambda(posL))]);
set(h,'FontSize',14);
subplot(3,1,3); h=plot(Xsol(:,posL+22)); 
set(h,'LineWidth',2);
hold on; h=plot(pos(1:S),x0(pos(1:S)),'ro'); 
set(h,'MarkerFaceColor','r');
set(gca,'FontSize',14);
h=text(50,0.65,['\lambda=',num2str(lambda(posL+22))]);
set(h,'FontSize',14);
% print -depsc2 Chapter_05_IRLS_Solutions.eps

figure(3); clf; 
h=semilogx(lambda,Xsol');
set(h,'LineWidth',2);
hold on;
h=semilogx(lambda,ones(101,1)*x0','--');
set(h,'LineWidth',2);
set(gca,'FontSize',14);
h=xlabel('\lambda'); 
set(h,'FontSize',14);
h=ylabel('Solution Entries');
set(h,'FontSize',14);
axis([min(lambda) max(lambda) -2 2]);
% print -depsc2 Chapter_05_IRLS_LambdaPath.eps

% Choosing the best lambbda and showing the IRLS there

lamOPT=lambda(posL);
Sol=zeros(m,101);
xout=zeros(m,1);
Sol(:,1)=xout;
XX=eye(m);
Err=zeros(101,2);
Err(1,1)=(xout-x0)'*(xout-x0); 
Err(1,2)=0.5*(A*xout-b)'*(A*xout-b)+lamOPT*sum(abs(xout)); 
for k=1:1:100
    xout=inv(2*lamOPT*XX+AA)*Ab;
    XX=diag(1./(abs(xout)+1e-10));
    Sol(:,k+1)=xout;
    Err(k+1,1)=(xout-x0)'*(xout-x0);
    Err(k+1,2)=0.5*(A*xout-b)'*(A*xout-b)+lamOPT*sum(abs(xout)); 
end;

figure(4); clf; 
semilogx(0:1:100,Sol');
hold on;
h=semilogx(0:1:100,ones(101,1)*x0',':'); 
set(h,'LineWidth',3);
h=xlabel('Iterations'); 
set(h,'FontSize',14);
h=ylabel('Solution Entries');
set(h,'FontSize',14);
set(gca,'FontSize',14);
% print -depsc2 Chapter_05_IRLS_Path.eps

figure(5); clf; 
h=semilogx(0:1:100,Err(:,2),'k'); 
set(h,'LineWidth',2);
h=xlabel('Iterations'); 
set(h,'FontSize',14);
h=ylabel('Function Value');
set(h,'FontSize',14);
set(gca,'FontSize',14);
grid on;
axis([1 100 1 2.2]);
% print -depsc2 Chapter_05_IRLS_Convergence.eps
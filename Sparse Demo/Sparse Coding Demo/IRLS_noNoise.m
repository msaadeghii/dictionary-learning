% Figure - NONE.  This code is from M. Elad's homepage.
% =========================================
% This script shows the IRLS (exact case) in practice, showing the 
% result after each iteration. Change s (line 7 ion the code) to see the 
% difference in behavior).
% Comment: Notice that this differs from the IRLS for the noisy case
%                 as described in Chapter 5.

% Prepare the problem Ax=b with sparse solution x
n=50; m=100; s=5; 
A=randn(n,m); 
A=A*diag(1./sqrt(diag(A'*A))); 
x_true=zeros(m,1); 
pos=randperm(m); 
x_true(pos(1:s))=randn(s,1); 
b=A*x_true;

% Solve the damn thing using IRLS
x=ones(m,1); % initialization
XX=eye(m); 

p=1; % refers to L0
p=0.5; % refers to L1

for k=1:1:20, 
    XX2=XX*XX; 
    x=XX2*A'*pinv(A*XX2*A')*b;
    XX=diag(abs(x).^p); 
    disp(['Iteration ',num2str(k),...
            ':    ||Ax-b||^2=',num2str((A*x-b)'*(A*x-b)),...
            ',    ||x-x_t||^2=',num2str((x-x_true)'*(x-x_true)),...
            '     ||x||_1=',num2str(sum(abs(x))),... 
            '     prod(x)=',num2str(exp(sum(log(abs(x)))/m))]);
    figure(1); clf; 
    pos=find(x_true ~=0); 
    h=plot(pos,x_true(pos),'or'); 
    set(h,'MarkerFaceColor','r'); 
    hold on; 
    stem(x,'b'); 
    pause  
end;
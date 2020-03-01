% This function reproduces the simulation results on synthetic dictionary
% recovery reported in the following paper:

%  Mostafa Sadeghi, Massoud Babaie-Zadeh, and Christian Jutten, 
%  “Dictionary learning for sparse representation: A novel approach,” 
%  IEEE Signal Processing Letters, vol. 20, no. 12, pp. 1195-1198, December 2013.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Mostafa Sadeghi
%     Electrical Engineering Department,
%     Sharif University of Technology,
%     Tehran, IRAN.
%     E-mail: m.saadeghii@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

K=50;                               % number of atoms
L=1500;                              % number of training signals
n=20;                             % signals dimension 
s=3;                                % number of non-zero coefficients

niter=30;                           % total number of runs
SNRdB=30;                           % noise level
param.numIteration=40;              % number of iterations in DL

ratio1=zeros(1,param.numIteration);
ratio2=zeros(1,param.numIteration);
ratio3=zeros(1,param.numIteration);

ratio4=zeros(1,param.numIteration);
ratio5=zeros(1,param.numIteration);
ratio6=zeros(1,param.numIteration);

tim1=0;tim2=0;tim3=0;
tim4=0;tim5=0;tim6=0;

for i=1:niter
[Do, Y, Xo] = generateSyntheticDictionaryAndData(L, s, n, K, SNRdB);

ind=randperm(L);
D=Y(:,ind(1:K));
D = normcols(D);

params.data = Y;
params.Tdata = s;
params.dictsize = K;
params.iternum = param.numIteration;
params.initdict=D;
params.trud=Do;
params.memusage = 'high';
params.innit=3;
params.dsp=1;
params.comperrdata=0;

% Dictionary learning methods

tic;[D1,g1,err1,ra1,rat1] = MOD_DL(params,'');t1=toc;tim1=tim1+t1;
ratio1=ratio1+rat1;

tic;[D2,g2,err2,ra2,rat2] = SGK_DL(params,'');t2=toc;tim2=tim2+t2;
ratio2=ratio2+rat2;

tic;[D3,g3,err3,ra3,rat3] = MDU_DL(params,'');t3=toc;tim3=tim3+t3;
ratio3=ratio3+rat3;

display(['Iteration: ',num2str(i),'/',num2str(niter), ' MOD-DL: ',num2str(rat1(end)),'... SGK-DL: ',num2str(rat2(end)),'... MDU-DL: ',num2str(rat3(end))])

% Dictionary Learning using the convex approximation approach

tic;[D4,g4,err4,ra4,rat4] = NewMOD(params,'');t4=toc;tim4=tim4+t4;
ratio4=ratio4+rat4;

tic;[D5,g5,err5,ra5,rat5] = NewSGK(params,'');t5=toc;tim5=tim5+t5;
ratio5=ratio5+rat5;

tic;[D6,g6,err6,ra6,rat6] = NewMDU(params,'');t6=toc;tim6=tim6+t6;
ratio6=ratio6+rat6;

display(['Iteration: ',num2str(i),'/',num2str(niter), ' NewMOD: ',num2str(rat4(end)),'... NewSGK: ',num2str(rat5(end)),'... NewMDU: ',num2str(rat6(end))])

end

tim1=tim1/niter;
tim2=tim2/niter;
tim3=tim3/niter;
tim4=tim4/niter;
tim5=tim5/niter;
tim6=tim6/niter;

rt1=ratio1/niter;
rt2=ratio2/niter;
rt3=ratio3/niter;
rt4=ratio4/niter;
rt5=ratio5/niter;
rt6=ratio6/niter;

display([' MOD: ',num2str(rt1(end)),'... SGK: ',num2str(rt2(end)),'... MDU: ',num2str(rt3(end))])
display('--------------Elapsed Times-------------')
display([' MOD: ',num2str(tim1),'... SGK: ',num2str(tim2),'... MDU: ',num2str(tim3)])

figure;plot(rt1,'k','linewidth',2);hold on;plot(rt2,'r','linewidth',2);plot(rt3,'b','linewidth',2);
legend('MOD','SGK','MDU');
grid on

display([' NewMOD: ',num2str(rt4(end)),'... NewSGK: ',num2str(rt5(end)),'... NewMDU: ',num2str(rt6(end))])
display('--------------Elapsed Times-------------')
display([' NewMOD: ',num2str(tim4),'... NewSGK: ',num2str(tim5),'... NewMDU: ',num2str(tim6)])

figure;plot(rt4,'k','linewidth',2);hold on;plot(rt5,'r','linewidth',2);plot(rt6,'b','linewidth',2);
legend('NewMOD','NewSGK','NewMDU');
grid on

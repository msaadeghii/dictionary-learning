% This function reproduces the simulation results on synthetic dictionary
% recovery reported in the following paper:

%  FatemeMostafa Sadeghi, Massoud Babaie-Zadeh, and Christian Jutten, 
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
s=6;                                % number of non-zero coefficients

niter=30;                           % total number of runs
SNRdB=100;                           % noise level
param.numIteration=200;              % number of iterations in DL

ratio1=zeros(1,param.numIteration);
ratio2=zeros(1,param.numIteration);
ratio3=zeros(1,param.numIteration);

tim1=0;tim2=0;tim3=0;

for i=1:niter
[Do, Y, Xo] = gererateSyntheticDictionaryAndData(L, s, n, K, SNRdB);

ind=randperm(L);
D=Y(:,ind(1:K));
D = normc(D);

params.data = Y;
params.Tdata = s;
params.dictsize = K;
params.iternum = param.numIteration;
params.initdict=D;
params.trud=Do;
params.memusage = 'high';
params.dsp=1;
params.comperrdata=0;

% Dictionary learning methods
 
%% KSVD

tic;[D1,g1,err1,ra1,rat1] = ksvd(params,'');t1=toc;tim1=tim1+t1;
ratio1=ratio1+rat1;

%% Proposed ADL

params.w=0.8;
params.initer=30;
params.tau=5e-3;
params.rho=1+1e-3;

params.lam=0.1;
params.w=0.85;

tic;[D2,g2,err2,ra2,rat2,iters2] = ADL(params,'');t2=toc;tim2=tim2+t2;
ratio2=ratio2+rat2;

%% MBHPA-DL

params.lam=0.1;
params.initer=1;
tic;[D3,g3,err3,ra3,rat3] = MBHPA(params,'');t3=toc;tim3=tim3+t3;
ratio3=ratio3+rat3;

%%

display(['Iteration: ',num2str(i),'/',num2str(niter), ' KSVD: '...
    ,num2str(rat1(end)),'... Proposed (ADL): ',num2str(rat2(end)),'  MBHPA-DL: ',num2str(rat3(end))])

end

tim1=tim1/niter;
tim2=tim2/niter;
tim3=tim3/niter;

rt1=ratio1/niter;
rt2=ratio2/niter;
rt3=ratio3/niter;

save(['Results_s_',num2str(s),'.mat'],'rt1','rt2','rt3','tim1','tim2','tim3');

figure;plot(rt1,'b','linewidth',2);hold on;plot(rt2,'r','linewidth',2);
plot(rt3,'m','linewidth',2);
legend('KSVD','Proposed (ADL)','MBHPA-DL')
grid on


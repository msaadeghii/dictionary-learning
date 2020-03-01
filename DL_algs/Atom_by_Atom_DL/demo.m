% A demo of how to use the functions described in this package. This script compare
% different DL algorithms, proposed in [1] and [2] in dictionary recovery from training data.

%   [1] M. Sadeghi, M. Babaie-Zadeh, and C. Jutten,
%         “Learning Overcomplete Dictionaries Based on Atom-by-Atom Updating,” 
%         IEEE Transactions on Signal Processing, vol. 62, no. 4, pp. 883-891, February 2014.
%
%   [2] M. Sadeghi, M. Babaie-Zadeh, and C. Jutten, 
%         “Learning over-complete dictionaries based on parallel atom-updating,” 
%         in Proceedings of 23rd IEEE International Workshop on Machine Learning for Signal Processing (MLSP 2013), London, 2013.

%%

clear all;close all;clc

%%

K=50;                                      % number of atoms
L=1500;                                  % number of training signals
n=20;                                      % signal (atom) length
s=5;                                         % number of non-zeros
SNRdB=100;                          % signal-to-noise ratio (dB)

niter=5;                                  % number of Monte-Carlo trials

params.innit=3;                     % number of inner-loop iterations
params.lam=1e-1;                 % regularization parameter of the l_1 norm
params.initx=eye(K,L);         % initial coefficient matrix
params.iternum=100;          % number of DL iterations
params.Tdata = s;                  % number of non-zeros
params.dictsize = K;              % number of atoms
params.memusage = 'high';

ratio1=zeros(1,params.iternum);ratio2=zeros(1,params.iternum);
ratio3=zeros(1,params.iternum);ratio4=zeros(1,params.iternum);

tim1=0;tim2=0;
tim3=0;tim4=0;

for i=1:niter % A more efficient way is to use "parfor", but it requires to remove "params" and send inputs to the functions one by one
    
    % synthetic training data generation
    [Do, Y, Xo] = generateSyntheticDictionaryAndData(L, s, n, K, SNRdB);
    % initial dictionary
    ind=randperm(L);
    D=Y(:,ind(1:K));
    D = normcols(D);

    params.data = Y;
    params.initdict=D;
    params.exact = 1;
    params.trud=Do;

    % Run DL algorithms
    
    tic;[D1,g1,err_ksvd1,ra1,rat1] = ksvd(params,'');t1=toc;tim1=tim1+t1;
    ratio1=ratio1+rat1;

    tic;[D2,g2,err_ksvd2,ra2,rat2] = PAU_DL(params,'');t2=toc;tim2=tim2+t2;
    ratio2=ratio2+rat2;
    
    tic;[D3,g3,err_ksvd3,ra3,rat3] = OS_DL(params,'');t3=toc;tim3=tim3+t3;
    ratio3=ratio3+rat3;
    
    tic;[D4,g4,err_ksvd4,ra4,rat4] = APrU_DL(params,'');t4=toc;tim4=tim4+t4;
    ratio4=ratio4+rat4;
    
    % Display the results
    
    display(['Iteration: ',num2str(i),'/',num2str(niter), ' ---- KSVD: ',num2str(rat1(end)),'... PAU-DL: ',num2str(rat2(end)),...
        '... OS-DL: ',num2str(rat3(end)),'... APrU_DL: ',num2str(rat2(end))])
    
end

tim1=tim1/niter;tim2=tim2/niter;
tim3=tim3/niter;tim4=tim4/niter;

rt1=ratio1/niter;rt2=ratio2/niter;
rt3=ratio3/niter;rt4=ratio4/niter;

%% Plot the results

itr1=1:params.iternum;
plot(itr1,rt1,'m',itr1,rt2,'r',itr1,rt3,'b',itr1,rt4,'c','linewidth',2.5);
set(gca,'fontsize',16);
h=legend('K-SVD','PAU-DL','OS-DL','APrU-DL');set(h,'fontsize',16);grid on
h1=xlabel('Iteration');set(h1,'FontSize',16);h1=ylabel('Successful recovery (%)');set(h1,'FontSize',16);
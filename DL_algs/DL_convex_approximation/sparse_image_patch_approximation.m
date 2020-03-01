% This function reproduces the simulation results on sparse image pathes approximation
% reported in the following paper:

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

Z = genpath(pwd);
addpath(Z);

sdata=load('data_image_patch.mat');
Dataa=sdata.Data;
Testdata=sdata.TestData;

% Set parameters
ddat=Dataa(:,1:20000);
imsize2=64;
params.data = ddat;
params.testdata=Testdata(:,1:5000);
params.Tdata = round(imsize2/10);
params.dictsize = 4*imsize2;
params.iternum = 50;
params.exact = 0;
params.memusage = 'high';
params.innit=3;
params.dsp=0;
params.comperrdata=1;

ind=randperm(size(ddat,2));
D=Dataa(:,ind(1:params.dictsize));
D = normcols(D);
params.initdict=D;


% Dictionary learning methods

tic;[D1,g1,err1,ra1,rat1] = MOD_DL(params,'i');t1=toc;
tic;[D2,g2,err2,ra2,rat2] = SGK_DL(params,'i');t2=toc;
tic;[D3,g3,err3,ra3,rat3] = MDU_DL(params,'i');t3=toc;

figure;plot(err1,'k','linewidth',2);hold on;plot(err2,'r','linewidth',2);plot(err3,'b','linewidth',2);
legend('MOD_DL','SGK_DL','MDU_DL');title('Training RMSE');
grid on

figure;plot(ra1,'k','linewidth',2);hold on;plot(ra2,'r','linewidth',2);plot(ra3,'b','linewidth',2);
legend('MOD-DL','SGK-DL','MDU-DL');title('Testing RMSE');
grid on

% Dictionary Learning using the convex approximation approach

tic;[D4,g4,err4,ra4,rat4] = NewMOD(params,'i');t4=toc;
tic;[D5,g5,err5,ra5,rat5] = NewSGK(params,'i');t5=toc;
tic;[D6,g6,err6,ra6,rat6] = NewMDU(params,'i');t6=toc;

figure;plot(err4,'k','linewidth',2);hold on;plot(err5,'r','linewidth',2);plot(err6,'b','linewidth',2);
legend('NewMOD','NewSGK','NewMDU');title('Training RMSE');
grid on

figure;plot(ra4,'k','linewidth',2);hold on;plot(ra5,'r','linewidth',2);plot(ra6,'b','linewidth',2);
legend('NewMOD','NewSGK','NewMDU');title('Testing RMSE');
grid on


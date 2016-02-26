% Reference:
%  Jianbo Yang Xuejun Liao, Xin Yuan, Patrick Llull, David J. Brady, Guillermo Sapiro and Lawrence Carin 
% ¡°Compressive Sensing by Learning a Gaussian Mixture Model from Measurements ¡± 
%  IEEE Transactions on Image Processing  (TIP)  vol. 24, no. 1, pp. 106-119, 2015.
% 
% Author:
%  Jianbo Yang,  the department of Electrical and Computer Engineering, Duke University
%  https://sites.google.com/site/jbysite/
% 
% Version history:
%  version 1.0 (July 10, 2014)
%  version 1.1 (Mar. 24, 2015)
 


clc
clear
close all
% --- Warning
% warning('off','all');

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);


filename = 'datasets/Kobe40_2_gray';  fname = 'Kobe40_2';
load(filename);
[Row,Col,nF] = size(Xtst);
Xtst = Xtst./(max(Xtst(:)))*255;
para.filename = filename;
para.fname = fname;

T = 8;  % How many frames collapsed to 1 measurement
M = 4;   % How many measurements to be used to reconstruct
para.patchSize = [4 4 T];
para.delta1 = 1;
para.delta2 = 1;
para.iteMax = 10; % How many iterations need to run for MMLE-GMM and MMLE-MFA. 
para.C = 8;

Xtst = Xtst(:,:,1:T*M);
C = binornd(1,0.5,[Row, Col,T]);
shift = 1;
for t=2:T
    C(:,1+(t-1)*shift:Col,t) = C(:,1+(t-2)*shift:Col-shift,t-1);
end
Y = zeros(Row,Col,M);
for t = 1:M
    Y(:,:,t) = sum(Xtst(:,:,(t-1)*T+(1:T)).*C,3);
end

model = [];
save('data_global','Y', 'C','Xtst');
para.flag_savemodel = true;
para.initType = 'gmmTrain'; % 1. psudoInverse; 2. loadFromFile; 3. preTrain
for nAlgorithm = 3
    tic
    switch nAlgorithm
        case 1 % Method GMM-TP
            para.method = 'GMM';                                            addpath('GMM/')
            [Xpre,model] = fGMM(para);                                      rmpath('GMM/')
        case 2 % Method MMLE-GMM
            para.method = 'GMM_FR_block_gmmTrain';                          addpath('GMM_FR/')
            [Xpre,model] = fGMM_FR(para);                                   rmpath('GMM_FR/')
        case 3 % Method MMLE-MFA
            para.method = 'GMM_LR_block_gmmTrain';                          addpath('GMM_LR/')
            [Xpre,model] = fGMM_LR(para);                                   rmpath('GMM_LR/')
    end
    time = toc;
    [psnr, ssim]= saveResults(Xpre, Row, Col, M, T,Y,Xtst,filename,para.method,time,model,nAlgorithm,para.C);
end


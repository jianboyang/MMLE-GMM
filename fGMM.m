function [Xpre,model] = fGMM(para)
load('data_global');
[para.Row, para.Col, para.T] = size(C);
para.M = size(Y,3);
patchSize = para.patchSize;
n = prod(patchSize);
delta1 = para.delta1;
delta2 = para.delta2;


para.R = 1e-6;
para.C = 20;



% A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);
% load(['GMM/' 'model_P' num2str(patchSize(1)) 'X' num2str(patchSize(2)) 'X' num2str(patchSize(3)) '_C' num2str(para.C)]);
% % load('GMM/model_P4X4X8_C20');
% % load('GMM/model_P4X4X14_C20');
% % load('GMM/model_P8X8X8_C20');
% model.Mu = Mu;
% model.Sig = Sig;
% model.pai = pai;
% for m = 1:para.M
%     y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);
%     X = GMM_3D_JBY(y, A, para, model);
%     Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.Row, para.Col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2);   
% end


blocksize = 64;
deltaBlock = blocksize/2;
Yorg = Y;
Corg = C;
Xtstorg = Xtst;
clear Y C Xtst;
Row = para.Row;
Col = para.Col;
T = para.T;
M = para.M;
para.Xprefile = 'xpreGMM';

% load(['GMM_model_C20_P4X4X8_20140331']);
load(['GMM/' 'model_P' num2str(patchSize(1)) 'X' num2str(patchSize(2)) 'X' num2str(patchSize(3)) '_C' num2str(para.C)]);
model.Mu = Mu;
model.Sig = Sig;
model.pai = pai;

Cb = video2patches_fast(Corg, blocksize, blocksize, para.T, deltaBlock, deltaBlock);
nBlock = size(Cb,2);
for m = 1:para.M
    Yb = image2patches_fast(Yorg(:,:,m), blocksize, blocksize, deltaBlock, deltaBlock);
    if ~isempty(Xtstorg)
        Xtstb = video2patches_fast(Xtstorg(:,:,(m-1)*para.T+1:m*para.T), blocksize, blocksize, para.T, deltaBlock, deltaBlock);
    end
   
    Xb = zeros(size(Cb,1), size(Cb,2));    
    fprintf('--------------------- m = %d ----------------------\n',m);
    C = reshape(Cb(:,1),blocksize,blocksize,T);
    [para.Row, para.Col, para.T] = size(C);
    para.n1 = para.patchSize(1);
    para.n2 = para.patchSize(2);
    para.delta1 = delta1;
    para.delta2 = delta2;
        
    matlabpool open local 4
    parfor ib = 1:nBlock      
        Y = reshape(Yb(:,ib),blocksize,blocksize);
        C = reshape(Cb(:,ib),blocksize,blocksize,T);
        if ~isempty(Xtstorg)
            Xtst = reshape(Xtstb(:,ib),blocksize,blocksize,T);
        else
            Xtst = [];
        end
        
        %     if mod(ib,10) == 1
        fprintf('%d of %d blocks:...\n',ib,nBlock);
        %     end
        
        A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);
        
        y = video2patches_fast(Y, para.patchSize(1), para.patchSize(2), 1, delta1, delta2);
        X = GMM_3D_JBY(y, A, para, model);
%         [X,model] = GMM_fly_fullrank(y, A, Xtst, para, [],xpreGMM);
        X = patches2video_fast(X, para.Row, para.Col, para.n1, para.n2, para.T, delta1, delta2);
        Xb(:,ib) = X(:);
    end
    matlabpool close
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(Xb, Row, Col, blocksize, blocksize, T, deltaBlock, deltaBlock);
end






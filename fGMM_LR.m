function [Xpre,model] = fGMM_LR(para)
load('data_global');
[para.Row, para.Col, para.T] = size(C);
para.M = size(Y,3);
patchSize = para.patchSize;
n = prod(patchSize);
delta1 = para.delta1;
delta2 = para.delta2;


blocksize = 64;
deltaBlock = blocksize/2;
Yorg = Y;
Corg = C;
Xtstorg = Xtst;
clear Y C Xtst;
Row = para.Row;
Col = para.Col;
para.burnin = 0;
para.rankS = 2;
para.Xprefile = [];
T = para.T;
M = para.M;


switch para.initType
    case 'gmmTrain'
        load([para.filename '_T' num2str(para.T) '_F' num2str(M) '_GMM']);
        XpreGMM = Xpre;   para.Xprefile = 'xpreGMM';
end

Cb = video2patches_fast(Corg, blocksize, blocksize, para.T, deltaBlock, deltaBlock);
nBlock = size(Cb,2);
for m = 1:para.M
    Yb = image2patches_fast(Yorg(:,:,m), blocksize, blocksize, deltaBlock, deltaBlock);
    if ~isempty(Xtstorg)
        Xtstb = video2patches_fast(Xtstorg(:,:,(m-1)*para.T+1:m*para.T), blocksize, blocksize, para.T, deltaBlock, deltaBlock);
    end
    switch para.initType
        case 'gmmTrain'
            Xpreb = video2patches_fast(XpreGMM(:,:,(m-1)*para.T+1:m*para.T), blocksize, blocksize, para.T, deltaBlock, deltaBlock);
    end
    
    Xb = zeros(size(Cb,1), size(Cb,2));
    fprintf('--------------------- m = %d ----------------------\n',m);
    C = reshape(Cb(:,1),blocksize,blocksize,T);
    [para.Row, para.Col, para.T] = size(C);
    para.n1 = para.patchSize(1);
    para.n2 = para.patchSize(2);
    para.delta1 = delta1;
    para.delta2 = delta2;
    
    for ib = 1:nBlock
        if strcmp(para.initType,'psudoInverse')
            xpreGMM{ib} = [];
        elseif strcmp(para.initType,'gmmTrain')
            xpreGMM{ib} = reshape(Xpreb(:,ib),blocksize,blocksize,T);
            xpreGMM{ib} = video2patches_fast(xpreGMM{ib}, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);
        else
            error('Wrong initType!');
        end
    end
    
    matlabpool open local 4
    parfor ib = 1:nBlock
        %             for ib = 1:nBlock
        Y = reshape(Yb(:,ib),blocksize,blocksize);
        C = reshape(Cb(:,ib),blocksize,blocksize,T);
        Xtst = reshape(Xtstb(:,ib),blocksize,blocksize,T);
        
        
        %     if mod(ib,10) == 1
        fprintf('%d of %d blocks:...\n',ib,nBlock);
        %     end
        
        A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);
        
        y = video2patches_fast(Y, para.patchSize(1), para.patchSize(2), 1, delta1, delta2);
        [X,model] = GMM_LR(y, A, Xtst, para, [],xpreGMM{ib});
        X = patches2video_fast(X, para.Row, para.Col, para.n1, para.n2, para.T, delta1, delta2);
        Xb(:,ib) = X(:);
    end
    matlabpool close
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(Xb, Row, Col, blocksize, blocksize, T, deltaBlock, deltaBlock);
end
model = [];

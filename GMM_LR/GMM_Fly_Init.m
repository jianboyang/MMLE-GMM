function [Mu, pai, F, subdim] = GMM_Fly_Init(para,model,CS,Y,xpreGMM)
[m,n] = size(CS{1});
delta1 = para.delta1;
delta2 = para.delta2;
rankS = para.rankS;
N = length(CS);
K = para.C;

subdim =  rankS*ones(1,K);
switch  para.initType
    case 'psudoInverse';
        for i=1:N
            [m,n] = size(CS{i});
            Xpre(:,i) = CS{i}'*inv(eye(m)*1e-5+CS{i}*CS{i}')*Y(:,i);
        end
    case 'gmmTrain'
        Xpre = xpreGMM;
    otherwise
        error('Choose the correct initialization');
end



if K == 1
    Mu = mean(Xpre,2);
    Sig{1} = cov(Xpre');
    pai = 1;
else
    if n < N
        options = statset('MaxIter',200);
        obj = gmdistribution.fit(Xpre',K,'Regularize',10^-4,'Options',options);
        
        pai = obj.PComponents;
        Mu = obj.mu';
        Sig = obj.Sigma;
        Sig = squeeze(num2cell(Sig,[1 2]));
    else
        [Mu,Sig,pai] = kmc(Xpre,K,'off');
        K = length(pai);
    end
end


for k=1:K
    [U,S,V] = svds(double(Sig{k}),subdim(k));
    S = sqrt(diag(S));
    F{k} = U*diag(S);
end

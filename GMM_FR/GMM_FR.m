function [Xpre, model] = GMM_FR(Y, CS, Xtst, para, model,xpreGMM)
stop_cond  = [1e-2,5e-6,para.iteMax];  % [1e-4,1e-4,100]
ite_N = stop_cond(3);
[m,n] = size(CS{1});
N = length(CS);
Xpre = zeros(n,N);
R = eye(m)*1e-5;
burnin = para.burnin;

t = 1;
converged = 0;
S = zeros(n,n);
llhi = zeros(N,1);
switch  para.initType
    case 'psudoInverse';
        for i=1:N
            Xpre(:,i) = CS{i}'*inv(eye(m)*1e-5+CS{i}*CS{i}')*Y(:,i);
        end
    case 'gmmTrain'
        Xpre = xpreGMM;
    otherwise
        error('Choose the correct initialization');
end



if para.C == 1    
    while ~converged && t < ite_N
        t = t+1;
        Mu = mean(Xpre,2);
        Sig = cov(Xpre');
        Sig = Sig + S/N;      
        
        S = zeros(n,n);
        for i = 1:N
            Q = Sig*CS{i}';
            P = inv(CS{i}*Q+R);
            P = (P+P')/2;
            Sig1 = Q*P;
            res = Y(:,i) - CS{i}*Mu;
            S = S + Sig - Sig1*Q';
            Xpre(:,i) = Mu + Sig1*res;
        end        
    end
    model.Mu = Mu;
    model.Sig = Sig;
else   % Soft Assignment
    Xpre = Xpre';
    label = k_means(Xpre, [], para.C);
    
    J = para.C;
    paii = full(sparse(1:N,label,1,N,J,N));  % Size of R is n X K
    S_bar = zeros(n,n,J);
    Sig = zeros(n,n,J);
    while ~converged && t < ite_N
        t = t+1;
        
        % Maximization
        nK = sum(paii);
        Mu = bsxfun(@times, paii'*Xpre, 1./nK');
        pai = sum(paii)/N;
        
        
        sqrtpaii = sqrt(paii);
        for j = 1:J
            Xo = bsxfun(@minus,Xpre,Mu(j,:));
            Xo = bsxfun(@times,Xo,sqrtpaii(:,j));
            Sig(:,:,j) = Xo'*Xo/nK(j) + eye(n)*(1e-6); % add a prior for numerical stability
            if t > burnin
                Sig(:,:,j) = Sig(:,:,j) + S_bar(:,:,j)/nK(j);
            end
        end
        Mu = Mu';
        
        
        % Expectation
        S_bar = zeros(n,n,J);
        for i = 1:N
            for j = 1:J
                mu1 = Mu(:,j);
                Q = Sig(:,:,j)*CS{i}';
                P = inv(CS{i}*Q+R);
                P = (P+P')/2;
                Sig1 = Q*P;
                res = Y(:,i) - CS{i}*mu1;
                tmp = -0.5*m*log(2*pi) + sum(log(diag(chol(P)))) - 0.5*sum(res.*(P*res),1);
                if t > burnin
                    S(:,:,j) = Sig(:,:,j) - Sig1*Q';
                end
                mu_M(:,j) = mu1 + Sig1*res;
                evd(j) = tmp + log(pai(j)+eps);
            end
            [paii(i,:),llhi(i)] = Design_Pai(evd);
            Xpre(i,:) = mu_M*paii(i,:)';
            if t > burnin
                for j = 1:J
                    S(:,:,j) = S(:,:,j)*paii(i,j);
                end
                S_bar = S_bar + S;
            end
        end
    end
    Xpre = Xpre';
    model.Mu = Mu;
    model.Sig = Sig;
    model.pai = pai;
end

function [paii,evm] = Design_Pai(evd)
c = size(evd,2);
mmax  = max(evd,[],2);
paii = exp(evd+repmat(-mmax,[1,c]))+eps;
evm = sum(paii,2);
paii = paii./repmat(evm,[1,c]);
evm = log(evm) + mmax;

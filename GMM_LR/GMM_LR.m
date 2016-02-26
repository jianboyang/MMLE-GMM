function [Xpre, model] = GMM_LR(Y, CS, Xtst,para, model,xpreGMM)
K = para.C;
stop_cond  = [1e-2,5e-3,para.iteMax];  % [1e-4,1e-4,100]
ite_N = stop_cond(3);
[m,n] = size(CS{1});
N = length(CS);
Xpre = zeros(n,N);
Gsigma = 1e-8;
Bsigma = 1e-8;
G = eye(n)*Gsigma;
B = eye(m)*Bsigma;



%% Initialization
[Mu, pai, F, subdim] = GMM_Fly_Init(para,model,CS,Y,xpreGMM);
if ite_N == 0
    model.pai = pai;
    model.Mu = Mu;
    model.F = F;
    model.G = G;
    model.B = B;
    return
end
llh = zeros(ite_N,1);


% Low rank GMM
for t=1:ite_N
    % compute p(k,z,x|z,Mu,F)
    mixp_new = zeros(1,K);
    meanXsum = zeros(n,K);
    meanZsum = cell(K,1);
    CorrXZsum = cell(K,1);
    CorrZZsum = cell(K,1);
    Psi_all = cell(K,1);
    for k=1:K
        meanZsum{k} = zeros(subdim(k),1);
        CorrXZsum{k} = zeros(n,subdim(k));
        CorrZZsum{k} = zeros(subdim(k),subdim(k));
        Z = cellfun(@(x) x*F{k},CS,'UniformOutput',false);
        Psi_all{k} = cat(3,Z{:});
    end
    llh(t) = 0;
    for i=1:N
        CorrZZ = cell(K,1);
        CorrXZ = cell(K,1);
        z = cell(K,1);
        log_rho = zeros(1,K);
        for k=1:K
            % using matrix inversion
            Delta = Bsigma + sum(CS{i}.^2,2)*Gsigma;
            %Psi = CS{i}*F(:,:,k);
            Psi = Psi_all{k}(:,:,i);
            invDeltaPsi = repmat(1./Delta,[1,subdim(k)]).*Psi;
            
            
            P = inv(eye(subdim(k))+Psi'*invDeltaPsi);   % (I +  Psi'*(Delta)^{-1}*Psi)^{-1}
            ydif   = Y(:,i) - CS{i}*Mu(:,k);
            
            term1 = P*invDeltaPsi';
            term2 = invDeltaPsi*term1;
            dist = ydif'*((1./Delta).*ydif) - ydif'*term2*ydif;
            [U,p] = chol(diag(1./Delta) - term2);
            logdetQ = 2*sum(log(diag(U)));
            log_rho(k) = log(pai(k)) - m*sum(log(2*pi))/2 -logdetQ/2 - dist/2;
            
                      
            z{k} = term1*ydif;
            x(:,k) = Mu(:,k) + F{k}*z{k} + Gsigma*CS{i}'*((1./Delta).*(ydif-Psi*z{k})) ;
                      
            CorrXZ{k} = (F{k} - Gsigma*CS{i}'*invDeltaPsi)*P + x(:,k)*z{k}';  % Adding x(:,k)*z{k}' is for computing muF
            CorrZZ{k} = P + z{k}*z{k}'; % Adding x(:,k)*z{k}' is for computing muF
        end
        base = max(log_rho);
        log_rho = log_rho - base;
        rho = exp(log_rho);
        rho = rho/sum(rho);
        
        Xpre(:,i) = x*rho';        
        mixp_new = mixp_new + rho;
        meanXsum = meanXsum + repmat(rho,[n,1]).*x;
        for k=1:K
            meanZsum{k}  = meanZsum{k} + repmat(rho(k),[subdim(k),1]).*z{k};
            CorrXZsum{k} = CorrXZsum{k} + rho(k)*CorrXZ{k};
            CorrZZsum{k} = CorrZZsum{k} + rho(k)*CorrZZ{k};
        end
        llh(t) = llh(t) + (log(sum(exp(log_rho))) + base)/N;
    end
    % update (Mu,F)
    pai = mixp_new/sum(mixp_new);
    for k=1:K,
        muF = [meanXsum(:,k),CorrXZsum{k}]*inv([mixp_new(k),meanZsum{k}'; meanZsum{k}, CorrZZsum{k}]);
        Mu(:,k) = muF(:,1);
        F{k} = muF(:,2:end);
    end
    %% check convergence
    if (t>3 & mean(abs(diff(llh(end-3:end))))/abs(llh(end))<stop_cond(2))...
            | 0%(t>1 & model.psnr(end)-model.psnr(end-1)<0)
        break;
    end
end

model.pai = pai;
model.Mu = Mu;
model.F = F;
model.llh = llh;
model.G = G;
model.B = B;




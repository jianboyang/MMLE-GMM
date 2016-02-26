function Xpre = GMM_3D_JBY(Y, CS, para,model)
T = para.T;
Mu = model.Mu;
Sig = model.Sig;
pai = model.pai;

[d,Nt] = size(Y);
Ri = para.R*eye(d);

%% Inversion
% GMM Inversion
Xpre = zeros(size(CS{1},2),Nt);

for i = 1:Nt
    [pai1,temp,~] = GMM_Inference_Jianbo(Y(:,i), CS{i}, Ri, Mu, Sig, pai);
    Xpre(:,i) = temp;
end



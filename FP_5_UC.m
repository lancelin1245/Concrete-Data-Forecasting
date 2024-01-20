%% UC Model
Vtau = 1; %% initial condition
omega2 = .04; %% fix omega^2
f = @(para) loglike_UC(para,omega2,y,Vtau);
parahat = fminsearch(f, [0.5, 0.5]);
phihat = parahat(2); sigma2hat = parahat(1);
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T)*phihat;
invOmega = sparse(1:T,1:T,[1/Vtau 1/omega2*ones(1,T-1)]);
HinvOmegaH = H'*invOmega*H;
K = speye(T)/sigma2hat + HinvOmegaH;
tauhat = K\(y/sigma2hat);

%% 
Vtau = 1; omega2 = 0.04; %% fix omega^2
h=1;
yhatUC = zeros(T-h-T0+1,1);
ytph = y(T0+h:end); % observed y_{t+h}
parahat = [1,1];
for t = T0:T-h
    yt = y(1:t);
    f = @(para) loglike_UC(para,omega2,yt,Vtau);
    parahat = fminsearch(f,parahat);
    phihat = parahat(2); sigma2hat = parahat(1);
    H = speye(t) - spdiags(ones(t-1,1),-1,t,t)*phihat;
    invOmega =sparse(1:t,1:t,[1/Vtau 1/omega2*ones(1,t-1)],t,t); 
    HinvOmegaH = H'*invOmega*H;
    K = speye(t)/sigma2hat + HinvOmegaH;
    tauhat = K\(yt/sigma2hat);
    yhatUC(t-T0+1) = tauhat(end); %store the forecasts
end
MSFE_6_UC = mean((ytph-yhatUC).^2);

%% UC loglikelihood
function ell = loglike_UC(para,omega,y,Vtau)
sigma = zeros(2,1);
sigma(1)= para(1); sigma(2)= para(2);
T = length(y);
%compute the MLE for tau given sig and omega
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T)*sigma(2);
invOmega = sparse(1:T,1:T,[1/Vtau 1/omega*ones(1,T-1)]); HinvOmegaH = H'*invOmega*H;
K = speye(T)/sigma(1) + HinvOmegaH;
tauhat = K\(y/sigma(1));
err = (y-tauhat)'*(y-tauhat)/sigma(1) + tauhat'*HinvOmegaH*tauhat;
ell = -T/2*log(sigma(1)) - (T-1)/2*log(omega) - .5*err;
ell = -ell;
end
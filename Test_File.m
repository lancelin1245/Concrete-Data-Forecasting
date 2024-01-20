%% Question 5 - IMA models
%% Question 5a, IMA(1,q) model
h = 1; % h-step-ahead forecast
q = 1; %the MA(q) model
yhatMA = zeros(T-h-T0+1,1); % forecasts
ytph = y(T0+h:end); % observed y_{t+h}
for t = T0:T-h
    yt = y(1:t); Lt = spdiags(ones(t-h+1,1),-1,t-h+1,t-h+1); It = speye(t);
    f = @(para) loglike_MAq(para,q,yt,y0); %para = (psi, sigma)
    if t == T0
        parahat = fminsearch(f, [zeros(1,q), 1]);
    else
        parahat = fminsearch(f, parahat);
    end
    psihat = parahat(1:q); sigmahat = parahat(q+1);
    Gamma = It;
    for k = 1:q
        Gamma = Gamma + psihat(k)*Lt^k;  %build Gamma depending on q
    end
    Xt = Lt*yt; Xt(1) = y(1) - y0(end);
    ehat = Gamma\(yt - Xt);
    yhatMA(t-T0+1,:) = y(t) + (flip(psihat))*ehat(end-q+1:end); %forecast
end
MSFE = mean((ytph-yhatMA).^2);
fprintf('Q5) IMA(1,q) model; q = %i, MSFE = %i\n', q, MSFE);
%fprintf('para: mu = %i, psi= %i', muhat, psihat);

%% Function for Question 5 - IMA(1,q) loglikelihood
function ell = loglike_MAq(para,q,y,y0) %para = (psi,sigma)
psi = para(1:q); sigma = para(q+1);
T = length(y); L = spdiags(ones(T,1),-1,T,T); I = speye(T); m =length(y0);
Gamma = I;
for k = 1:q
    Gamma = Gamma + psi(k)*L^k;  %build Gamma depending on q
end
X = L*y; X(1) = y(1) - y0(end);
ell = -(T/2)*log(2*pi*(sigma^2)) - ((((y - X)')/(Gamma'))*(Gamma\(y - X)))/(2*sigma^2);
ell = -ell;
end

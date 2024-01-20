%% ARX(1) with AR(2) errors

L = spdiags(ones(T,1),-1,T,T); I = speye(T);
X = [ones(T,1), (1:T)', L*y]; 

ytph_ARX1_AR2errors = y(T0+h1:end); % observed y {t+h} 
yhat_ARX1_AR2errors = zeros(T-h1-T0+1,1); % ARX(1) forecast
for t = T0:T-h1
	yt = y(h1:t); Lt = spdiags(ones(t-h1+1,1),-1,t-h1+1,t-h1+1);
    f = @(psi) cloglike_ARX1AR2(psi,yt,y0);
    if t == T0
        [psi, ~] = fminsearch(f,[-0.2, -0.2]); %first go, use initial value psi = (-0.2, -0.2)
    else
        [psi, ~] = fminsearch(f,psi); %afterwards, use previous psi as inital
    end
    Xt = [ones(t-h1+1,1), (h1:t)', Lt*yt]; 
    Xt(1,3) = y0(end);  %include lags   
    H = I - psi(1)*L - psi(2)*(L^2);
    betahat = (X'*(H'*H)*X)\(X'*(H'*H)*y);
    ehat = yt - Xt*betahat; %estimate e
    yhat_ARX1_AR2errors(t-T0+1,:) = [1, t+1, y(t)]*betahat + ehat(end)*psi(1) + ehat(end-1)*psi(2); %store the forecasts
end
ytph = y(T0+h1:end); % observed y_{t+h}
MSFE_3_4_ARX1_AR2errors = mean((ytph_ARX1_AR2errors-yhat_ARX1_AR2errors).^2);


%% Loglikelihood function for ARX(1) model with AR(2) errors
function ell = cloglike_ARX1AR2(para,y,y0)
psi = para;
T = length(y); L = spdiags(ones(T,1),-1,T,T); I = speye(T); m =length(y0);
X = [ones(T,1), (1:T)', L*y]; X(1,3) = y0(end);   %include lags
H = I - psi(1)*L - psi(2)*(L^2);
betahat = (X'*H'*H*X)\(X'*H'*H*y);
sigma2hat = (y - X*betahat)'*H'*H*(y - X*betahat)/T;
ell = -(T/2)*log(2*pi*sigma2hat) - (1/(2*sigma2hat))*(y - X*betahat)'*H'*H*(y - X*betahat);
ell = -ell;
end
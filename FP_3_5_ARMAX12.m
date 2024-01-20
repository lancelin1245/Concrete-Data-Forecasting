%% ARMAX(1,2)
ytph_ARMAX12 = y(T0+h1:end); % observed y {t+h} 
yhat_ARMAX12 = zeros(T-h1-T0+1,1); % ARX(1) forecast

for t = T0:T-h
	yt = y(h:t); Lt = spdiags(ones(t-h+1,1),-1,t-h+1,t-h+1); It = speye(t-h+1);
    f = @(psi) cloglike_ARX1MA2(psi,yt,y0);
    if t == T0
        [psi, ~] = fminsearch(f,[-0.1,-0.1]); %first go, use initial value psi = 0
    else
        [psi, ~] = fminsearch(f,[-0.1,-0.1]); %normally use old psi as the initial value, but this causes cascading errors here
    end
	Xt = [ones(t-h+1,1), (h:t)', Lt*yt]; 
    Xt(1,3) = y0(end); %include lags
	Gamma = It + psi(1)*Lt + psi(2)*(Lt^2);
    Gamma2 = Gamma*Gamma';
    betahat = (Xt'*(Gamma2\Xt))\(Xt'*(Gamma2\yt));
    uhat = Gamma\(yt - Xt*betahat); %estimate u
    yhat_ARMAX12(t-T0+1,:) = [1, t+1, y(t)]*betahat + uhat(end)*psi(1) + uhat(end-1)*psi(2); %store the forecasts
end
ytph = y(T0+h:end); % observed y_{t+h}
MSFE_ARMAX12 = mean((ytph_ARMAX12-yhat_ARMAX12).^2);

%% Loglikelihood function for ARMAX(1,2)
function ell = cloglike_ARX1MA2(para,y,y0)
psi = para;
T = length(y); L = spdiags(ones(T,1),-1,T,T); I = speye(T); m =length(y0);
X = [ones(T,1), (1:T)', L*y]; 
X(1,3) = y0(end);   %include lags
Gamma = I + psi(1)*L + psi(2)*(L^2);
Gamma2 = Gamma*Gamma';
prebetahat = (X'*(Gamma2\X))\(X'/Gamma2);
betahat = prebetahat*y;
%betahat = (X'*(Gamma2\X))\(X'*(Gamma2\y));
sigma2hat = (y - X*betahat)'*(Gamma2\(y - X*betahat))*(1/T);
ell = -(T/2)*log(2*pi*sigma2hat) - (1/(2*sigma2hat))*(y - X*betahat)'*(Gamma2\(y - X*betahat));
ell = -ell;
end
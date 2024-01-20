%% IMA(1,1), 1-step forecast
dely = y_d1;

yhat_MA1_4step = zeros(T-h1-T0+1,1); %% IMA(1,1) forecasts 
ima11_4step = @(psi) loglike_MA1(psi,dely(1:T0)); 
psihat_ima11_4step = fminsearch(ima11_4step,[0.5,0.5]);

for t = T0:T
    delyt = dely(1:t); 
    % find the MLE 
    ima11_4step = @(psi) loglike_MA1(psi,delyt); 
    psihat_ima11_4step = fminsearch(ima11_4step, psihat_ima11_4step); 
    % prediction 
    Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima11_4step(2); 
    uhat = Gamma\delyt; 
    % store the forecasts
    yhat_MA1_4step(t-T0+1,:) = y(t) + psihat_ima11_4step(2)*uhat(end);
end 

%% negative of the log lilkelihood for MA(1) 
% input: x = [sig, theta1]; y = data 
function ell = loglike_MA1(x,y) 
    sig = x(1); theta1 = x(2); 
    T = length(y); 
    A = speye(T); 
    B = spdiags(ones(T-1,1),[-1],T,T); 
    Gam = A + B*theta1; 
    Gam2 = Gam*Gam'; 
    ell = -T/2*log(2*pi*sig)-.5*log(det(Gam2))-.5/sig*y'*(Gam2\y); 
    ell = -ell;
end
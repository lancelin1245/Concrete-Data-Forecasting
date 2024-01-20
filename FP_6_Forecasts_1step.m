%% IMA(1,1), 1-step forecast
h1 = 1;
dely_1step = totaly(m+h1:end) - totaly(m:end-h1);

yhat_MA1_1step = zeros(T-h1-T0+1,1); %% IMA(1,1) forecasts 
ima11_1step = @(psi) loglike_MA1(psi,dely_1step(1:T0)); 
psihat_ima11_1step = fminsearchbnd(ima11_1step,[0.5,0.5]);

for t = T0:T
    delyt = dely_1step(1:t); 
    % find the MLE 
    ima11_1step = @(psi) loglike_MA1(psi,delyt); 
    psihat_ima11_1step = fminsearch(ima11_1step, psihat_ima11_1step); 
    % prediction 
    Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima11_1step(2); 
    uhat = Gamma\delyt; 
    % store the forecasts
    yhat_MA1_1step(t-T0+1,:) = y(t) + psihat_ima11_1step(2)*uhat(end);
end 

% yhat_t+1 = 1.1772
% Interval forecast = (1177200 -1.96*sqrt(psihat_ima11_1step(1) +
% psihat_ima11_1step(1)*psihat_ima11_1step(2)^2), 1177200 +
% 1.96*sqrt(psihat_ima11_1step(1) + psihat_ima11_1step(1)*psihat_ima11_1step(2)^2))
yhat_h2 = yhat_MA1_1step(end) + psihat_ima11_1step(2)*uhat;
yhat_h3 = yhat_h2(end) + psihat_ima11_1step(2)*uhat;
yhat_h4 = yhat_h3(end) + psihat_ima11_1step(2)*uhat;

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
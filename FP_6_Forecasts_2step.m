%% IMA(1,1), 2-step forecast
h2 = 2;
dely_2step = totaly(m+h2:end) - totaly(m:end-h2);

yhat_MA1_2step = zeros(T-h2-T0+1,1); %% IMA(1,1) forecasts 
ima11_2step = @(psi) loglike_MA1(psi,dely_2step(1:T0)); 
psihat_ima11_2step = fminsearch(ima11_2step,[0.5,0.5]);
% 
% for t = T0:T-h2
%     delyt = dely_2step(1:t); 
%     % find the MLE 
%     ima11_2step = @(psi) loglike_MA1(psi,delyt); 
%     psihat_ima11_2step = fminsearch(ima11_2step, psihat_ima11_2step); 
%     % prediction 
%     Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima11_2step(2); 
%     uhat = Gamma\delyt; 
%     % store the forecasts
%     yhat_MA1_2step(t-T0+1,:) = y(t) + psihat_ima11_2step(2)*uhat(end);
% end 

Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima11_1step(2); 
uhat = Gamma\delyt;
yhat_h2 = yhat_MA1_1step(end) + psihat_ima11_2step(2)*uhat;

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
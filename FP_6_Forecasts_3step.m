%% IMA(1,1), 1-step forecast
h3 = 3;
dely_3step = totaly(m+h3:end) - totaly(m:end-h3);

yhat_MA1_3step = zeros(T-h3-T0+1,1); %% IMA(1,1) forecasts 
ima11_3step = @(psi) loglike_MA1(psi,dely_3step(1:T0)); 
psihat_ima11_3step = fminsearch(ima11_3step,[0.5,0.5]);

% for t = T0:T
%     delyt = dely(1:t); 
%     % find the MLE 
%     ima11_3step = @(psi) loglike_MA1(psi,delyt); 
%     psihat_ima11_3step = fminsearch(ima11_3step, psihat_ima11_3step); 
%     % prediction 
%     Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima11_3step(2); 
%     uhat = Gamma\delyt; 
%     % store the forecasts
%     yhat_MA1_3step(t-T0+1,:) = y(t) + psihat_ima11_3step(2)*uhat(end);
% end 

Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima11_1step(2); 
uhat = Gamma\delyt;
yhat_h3 = yhat_MA1_2step(end) + psihat_ima11_3step(2)*uhat;

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


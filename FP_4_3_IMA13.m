%% IMA(1,3), 1-step forecast
yhat_MA3 = zeros(T-h1-T0+1,1); %% IMA(1,1) forecasts 
ytph = y(T0+h1:end); % observed y_{t+h} 
ima13 = @(psi) loglike_MA3(psi,dely(1:T0)); 
psihat_ima13 = fminsearch(ima13,[0.5,0.5,0.5,0.5]);

for t = T0:T-h1
    delyt = dely(1:t); 
    % find the MLE 
    ima13 = @(psi) loglike_MA3(psi,delyt); 
    psihat_ima13 = fminsearch(ima13, psihat_ima13); 
    % prediction 
    Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima13(2) + spdiags(ones(t-2,1),[-2],t,t)*psihat_ima13(3) + spdiags(ones(t-3,1),[-3],t,t)*psihat_ima13(4); 
    uhat = Gamma\delyt; 
    % store the forecasts
    yhat_MA3(t-T0+1,:) = y(t) + psihat_ima13(2)*uhat(end) + psihat_ima13(3)*uhat(end-1) + psihat_ima13(4)*uhat(end-2);
end 
MSFE_MA3 = mean((ytph-yhat_MA3).^2); % 0.4090

%% negative of the log lilkelihood for MA(3) 
% input: x = [sig, theta1, theta2, theta3]; y = data 
function ell = loglike_MA3(x,y) 
    sig = x(1); theta1 = x(2); theta2 = x(3); theta3 = x(4);
    T = length(y); 
    A = speye(T); 
    B = spdiags(ones(T-1,1),[-1],T,T); 
    C = spdiags(ones(T-2,1),[-2],T,T);
    D = spdiags(ones(T-3,1),[-3],T,T);
    Gam = A + B*theta1 + C*theta2 + D*theta3; 
    Gam2 = Gam*Gam'; 
    ell = -T/2*log(2*pi*sig)-.5*log(det(Gam2))-.5/sig*y'*(Gam2\y); 
    ell = -ell;
end
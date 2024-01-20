%% IMA(1,2), 1-step forecast
yhat_MA2 = zeros(T-h1-T0+1,1); %% IMA(1,1) forecasts 
ytph = y(T0+h1:end); % observed y_{t+h} 
ima12 = @(psi) loglike_MA2(psi,dely(1:T0)); 
psihat_ima12 = fminsearch(ima12,[0.5,0.5,0.5]);

for t = T0:T-h1
    delyt = dely(1:t); 
    % find the MLE 
    ima12 = @(psi) loglike_MA2(psi,delyt); 
    psihat_ima12 = fminsearch(ima12, psihat_ima12); 
    % prediction 
    Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima12(2) + spdiags(ones(t-2,1),[-2],t,t)*psihat_ima12(3); 
    uhat = Gamma\delyt; 
    % store the forecasts
    yhat_MA2(t-T0+1,:) = y(t) + psihat_ima12(2)*uhat(end) + psihat_ima12(3)*uhat(end-1);
end 
MSFE_MA2 = mean((ytph-yhat_MA2).^2); % 0.3998

%% IMA(1,2) 1-step plot
hold on
plot(y)
plot(T0 + 1:T, yhat_MA2)
legend('Actual Data','1-Step Ahead IMA(1,2) Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for IMA(1,2) Model', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)

%% negative of the log lilkelihood for MA(2) 
% input: x = [sig, theta1, theta2]; y = data 
function ell = loglike_MA2(x,y) 
    sig = x(1); theta1 = x(2); theta2 = x(3); 
    T = length(y); 
    A = speye(T); 
    B = spdiags(ones(T-1,1),[-1],T,T); 
    C = spdiags(ones(T-2,1),[-2],T,T);
    Gam = A + B*theta1 + C*theta2; 
    Gam2 = Gam*Gam'; 
    ell = -T/2*log(2*pi*sig)-.5*log(det(Gam2))-.5/sig*y'*(Gam2\y); 
    ell = -ell;
end
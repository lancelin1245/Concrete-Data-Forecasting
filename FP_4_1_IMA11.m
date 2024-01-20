%% IMA(1,1), 1-step forecast
dely = y_d1;

yhat_MA1 = zeros(T-h1-T0+1,1); %% IMA(1,1) forecasts 
ytph = y(T0+h1:end); % observed y_{t+h} 
ima11 = @(psi) loglike_MA1(psi,dely(1:T0)); 
psihat_ima11 = fminsearch(ima11,[0.5,0.5]);

for t = T0:T-h1
    delyt = dely(1:t); 
    % find the MLE 
    ima11 = @(psi) loglike_MA1(psi,delyt); 
    psihat_ima11 = fminsearch(ima11, psihat_ima11); 
    % prediction 
    Gamma = speye(t) + spdiags(ones(t-1,1),[-1],t,t)*psihat_ima11(2); 
    uhat = Gamma\delyt; 
    % store the forecasts
    yhat_MA1(t-T0+1,:) = y(t) + psihat_ima11(2)*uhat(end);
end 
MSFE_MA1 = mean((ytph-yhat_MA1).^2); % 5.0858e+09

%% IMA(1,1) 1-step plot
hold on
plot(y)
plot(T0 + 1:T, yhat_MA1)
legend('Actual Data','1-Step Ahead IMA(1,1) Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for IMA(1,1) Model', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)

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
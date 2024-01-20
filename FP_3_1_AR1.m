%% AR(1)
% 1-step
ytph_AR1_h1 = y(T0+h1:end); % observed y {t+h} 
yhat_AR1_h1 = zeros(T-h1-T0+1,1); % AR(1) forecast
for t = T0:T-h1
    yt = y(h1:t); 
    Xt = [ones(t-h1+1,1) y_lag1(1:t-h1+1)]; 
    betahat = (Xt'*Xt)\(Xt'*yt); 
    % store the forecasts 
    yhat_AR1_h1(t-T0+1,:) = [1 y(t)]*betahat;
end 
MSFE_AR1_h1 = mean((ytph_AR1_h1-yhat_AR1_h1).^2); % 6.4613e+09

%% AR(1) 1-step plot
hold on
plot(y)
plot(T0 + 1:T, yhat_AR1_h1)
legend('Actual Data','1-Step Ahead AR(1) Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for AR(1) Model', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)




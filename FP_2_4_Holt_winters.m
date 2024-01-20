%% Holt Winters 1-step
syhat_2_4_1step = zeros(T-h1-T0+1,1);
ytph_2_4_1step = y(T0+h1:end); % observed y_{t+h}
alpha = .3; beta = .3; gamma = .3; % smoothing parameters
s = 4; % periodicity of seasonality
Lt = mean(y(1:s)); bt = 0; St(1:s) = y(1:s)-Lt;
for t = s+1:T-h1
	newLt = alpha*(y(t)-St(t-s)) + (1-alpha)*(Lt+bt); % level
	newbt = beta*(newLt-Lt) + (1-beta)*bt; % slope
	St(t) = gamma*(y(t)-newLt) + (1-gamma)*St(t-s); % Seasonality
    yhat = newLt + h1*newbt + St(t+h1-s);
	Lt = newLt; bt = newbt; % update Lt and bt
	if t>= T0 % store the forecasts for t >= T0
		syhat_2_4_1step(t-T0+1,:) = yhat;
	end
end
MSFE_2_4_1step = mean((ytph_2_4_1step-syhat_2_4_1step).^2);

%% Holt-Winter 1-step plot
hold on
plot(1:T, y)
plot(T0 + 1:T, syhat_2_4_1step)
legend('Actual Data','1-Step Ahead Holt-Winter Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for Holt-Winter 1-step Forecast', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)
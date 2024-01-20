%% Random walk
% Random walk model without drift
% y_that = y_t-1 

y_rw = [];

for i = T0:T
    y_rw(length(y_rw) + 1) = y_lag1(i);
end

y_rw = y_rw';

MSFE_rw = mean((y(T0:end) - y_rw(1:end)).^2);

y_rw_forecast = y(end);

base = MSFE_rw;
%%
hold on
plot(y)
plot(T0:T, y_rw, '--')
legend('Actual Data','Random Walk Forecast', 'FontSize',25, 'Location','northwest')
hold off
title('Plot of random-walk forecast', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25) 

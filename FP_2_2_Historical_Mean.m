%% Historical mean 

y_hm = [];

for i = T0:T
    y_hm(length(y_hm) + 1) = mean(y(1:i));
end

y_hm = y_hm';

MSFE_hm = mean((y(T0:end) - y_hm(1:end)).^2);

y_hm_forecast = mean(y(1:T));

%%
hold on
plot(y)
plot(T0:T, y_hm, '--')
legend('Actual Data', 'Historical Mean Forecast', 'FontSize',25, 'Location','northwest')
hold off
title('Plot of historical-mean forecast', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25) 
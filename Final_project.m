%% Set up data
Concrete = readmatrix('Concrete_data.csv');
totaly = Concrete(1:118, 2);
y = Concrete(9:118, 2); % T0 = 8, 1995Q1 to 2022Q2
T = length(y); t = (1:T)'; total_T = length(totaly);
m=8;
y0 = totaly(1:m, 1); % 8 lags that can be used
y_lag1 = totaly(m:end-1, 1);
y_lag2 = totaly(m-1:end-2, 1);
y_lag3 = totaly(m-2:end-3, 1);
y_lag4 = totaly(m-3:end-4, 1);
y_lag5 = totaly(m-4:end-5, 1);
y_lag6 = totaly(m-5:end-6, 1);
y_lag7 = totaly(m-6:end-7, 1);
y_lag8 = totaly(1:end-8, 1);

h1 = 1;
T0 = 8;
% [h,pValue,stat,cValue] = adftest(y)
% [h,pValue,stat,cValue] = adftest(y_d1)
% base = MSFE_rw;
%% Plotting data
plot(y)
title('Plot of ready-mixed concrete from 1995Q1 to 2022Q2', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)

%% Stem plot
[cov_y,lags_y] = xcov(y, 50, 'coeff');
stem(lags_y,cov_y)
title('Plot of sample autocovariance for concrete data', 'FontSize',30)

%%
% First differenced
y_d1 = zeros(T-1, 1);

for t = 1:T
    if t == 1
    y_d1(t) = y(t) - y_lag1(1);
    else
    y_d1(t) = y(t) - y(t-1);
    end
end
%% 
plot(y_d1)
title('Plot of first-differenced concrete data', 'FontSize',30)

%%
[cov_y_d1,lags_y_d1] = xcov(y_d1, 50, 'coeff');
stem(lags_y_d1,cov_y_d1)   
title('Plot of sample autocovariance for first-differenced concrete data', 'FontSize',30)

%%
plot(log(y))
title('Plot of logged concrete data')

%%
log_y = log(y);
log_y_d1 = zeros(T-1, 1);

for t = 1:T
    if t == 1
    log_y_d1(t) = log_y(t) - log(y_lag1(1));
    else
    log_y_d1(t) = log_y(t) - log_y(t-1);
    end
end
plot(log_y_d1)
title('Plot of first-differenced logged Concrete data', 'FontSize',30)

%%
[cov_log_y_d1,lags_log_y_d1] = xcov(log_y_d1, 50, 'coeff');
stem(lags_log_y_d1,cov_log_y_d1) 
title('Plot of sample autocovariance for first-differenced logged concrete data', 'FontSize',30)

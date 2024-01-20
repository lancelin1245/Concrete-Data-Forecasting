%% Seasonal dummies
Q = Concrete(9:118, 5);

% Dummy variables for quarters
D1 = (Q == 1); D2 = (Q == 2); D3 = (Q == 3); D4 = (Q == 4);
%% 1st specification : yt = a0 + a1*t + a1*D1 + e
% 1 step
syhat_1st_1step = zeros(T-h1-T0+1,1);
ytph_1st_1step = y(T0+h1:end); % observed y_{t+h}
for t = T0:T-h1
    yt = y(1:t);
    D1t = D1(1:t); D2t = D2(1:t);
    D3t = D3(1:t); D4t = D4(1:t);
    Xt = [ones(t,1) (1:t)' D1t];
    beta1 = (Xt'*Xt)\(Xt'*yt);
    yhat1 = [1 t+h1 D1(t+h1)]*beta1;
    syhat_1st_1step(t-T0+1) = yhat1;
end
MSFE_1st_2_3 = mean((ytph_1st_1step - syhat_1st_1step).^2);

%%
hold on
plot(1:T, y)
plot(T0 + 1:T, syhat_1st_1step, '--')
legend('Actual Data','1-Step Ahead Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for First Specification', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)

%% 2nd specification : yt = a0 + a1*t + a2*D2 + e
% 1 step
syhat_2nd_1step = zeros(T-h1-T0+1,1);
ytph_2nd_1step = y(T0+h1:end); % observed y_{t+h}
for t = T0:T-h1
    yt = y(1:t);
    D1t = D1(1:t); D2t = D2(1:t);
    D3t = D3(1:t); D4t = D4(1:t);
    Xt = [ones(t,1) (1:t)' D2t];
    beta1 = (Xt'*Xt)\(Xt'*yt);
    yhat1 = [1 t+h1 D2(t+h1)]*beta1;
    syhat_2nd_1step(t-T0+1) = yhat1;
end
MSFE_2nd_2_3 = mean((ytph_2nd_1step - syhat_2nd_1step).^2);

%%
hold on
plot(1:T, y)
plot(T0 + 1:T, syhat_2nd_1step)
legend('Actual Data','1-Step Ahead Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for Second Specification', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)
%% 3rd specification : yt = a0 + a1*t + a3*D3 + e
% 1 step
syhat_3rd_1step = zeros(T-h1-T0+1,1);
ytph_3rd_1step = y(T0+h1:end); % observed y_{t+h}
for t = T0:T-h1
    yt = y(1:t);
    D1t = D1(1:t); D2t = D2(1:t);
    D3t = D3(1:t); D4t = D4(1:t);
    Xt = [ones(t,1) (1:t)' D3t];
    beta1 = (Xt'*Xt)\(Xt'*yt);
    yhat1 = [1 t+h1 D3(t+h1)]*beta1;
    syhat_3rd_1step(t-T0+1) = yhat1;
end
MSFE_3rd_2_3 = mean((ytph_3rd_1step - syhat_3rd_1step).^2);

%%
hold on
plot(1:T, y)
plot(T0 + 1:T, syhat_3rd_1step)
legend('Actual Data','1-Step Ahead Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for Third Specification', 'FontSize',30)
    ylabel('Volume (m^{3})', 'FontSize',25)

%% 4th specification : yt = a0 + a1*t + a4*D4 + e
% 1 step
syhat_4th_1step = zeros(T-h1-T0+1,1);
ytph_4th_1step = y(T0+h1:end); % observed y_{t+h}
for t = T0:T-h1
    yt = y(1:t);
    D1t = D1(1:t); D2t = D2(1:t);
    D3t = D3(1:t); D4t = D4(1:t);
    Xt = [ones(t,1) (1:t)' D4t];
    beta1 = (Xt'*Xt)\(Xt'*yt);
    yhat1 = [1 t+h1 D4(t+h1)]*beta1;
    syhat_4th_1step(t-T0+1) = yhat1;
end
MSFE_4th_2_3 = mean((ytph_4th_1step - syhat_4th_1step).^2);

%%
hold on
plot(1:T, y)
plot(T0 + 1:T, syhat_4th_1step)
legend('Actual Data','1-Step Ahead Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for Fourth Specification', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)

%% 5th specification : yt = a1*t + a1*D1 + a2*D2 + a3*D3+ a4*D4 + e
% 1 step
syhat_5th_1step = zeros(T-h1-T0+1,1);
ytph_5th_1step = y(T0+h1:end); % observed y_{t+h}
for t = T0:T-h1
    yt = y(1:t);
    D1t = D1(1:t); D2t = D2(1:t);
    D3t = D3(1:t); D4t = D4(1:t);
    Xt = [(1:t)' D1t D2t D3t D4t];
    beta1 = (Xt'*Xt)\(Xt'*yt);
    yhat1 = [t+h1 D1(t+h1) D2(t+h1) D3(t+h1) D4(t+h1)]*beta1;
    syhat_5th_1step(t-T0+1) = yhat1;
end
MSFE_5th_2_3 = mean((ytph_5th_1step - syhat_5th_1step).^2);

%%
hold on
plot(1:T, y)
plot(T0 + 1:T, syhat_5th_1step)
legend('Actual Data','1-Step Ahead Forecast', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for Fifth Specification', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)
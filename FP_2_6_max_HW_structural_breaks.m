%% Structural break for GFC
% GFC = Concrete(:, 6);
% 
v = [0,0,0];
optimisers_2_6_GFC_1step = fminsearchbnd(@holtoptimiser_2_6_GFC_1step, v, [0, 0, 0], [1, 1, 1]);
optimisers_2_6_GFC_1step

%%
% 1-step
% alpha = 0.4391 beta = 0.1636 gamma = 0.2123 MSFE = 4.1690e+09

MSFE_2_6_GFC_1step = holtoptimiser_2_6_GFC_1step([0.4391, 0.1636, 0.2123])

%%
syhat_2_6_GFC_1step = zeros(T-h1-T0+1,1);
ytph_2_6_GFC_1step = y(T0+h1:end); % observed y_{t+h}
alpha = 0.4391; beta = 0.1636; gamma = 0.2123; % smoothing parameters
s = 4; % periodicity of seasonality
Lt = mean(y(1:s)); bt = 0; St(1:s) = y(1:s)-Lt;
GFC = 56; % start of GFC period 2008Q4
psi = .7; % parameter that exponentially becomes smaller as period after GFC increases
for t = s+1:T-h1
    if (t < GFC)
	    newLt = alpha*(y(t)-St(t-s)) + (1-alpha)*(Lt+bt); % level
	    newbt = beta*(newLt-Lt) + (1-beta)*bt; % slope
        St(t) = gamma*(y(t)-newLt) + (1-gamma)*St(t-s); % Seasonality
        yhat = newLt + h1*newbt + St(t+h1-s);
	    Lt = newLt; bt = newbt; % update Lt and bt
    elseif (t >= GFC)
        newLt = alpha*(y(t)-St(t-s)) + (1-alpha)*(Lt+bt); % level
	    newbt = beta*(newLt-Lt) + (1-beta)*bt; % slope
        St(t) = gamma*(y(t)-newLt) + (1-gamma)*St(t-s); % Seasonality
        yhat = newLt + h1*newbt + St(t+h1-s) - (15595)*(psi^(t+1-GFC));
	    Lt = newLt; bt = newbt; % update Lt and bt
    end
    
    if t>= T0 % store the forecasts for t >= T0
		syhat_2_6_GFC_1step(t-T0+1,:) = yhat;
 	end
end

MSFE_2_6_GFC_1step = mean((ytph_2_6_GFC_1step-syhat_2_6_GFC_1step).^2);

%% Holt-Winter 1-step plot
hold on
plot(1:T, y)
plot(T0 + 1:T, syhat_2_6_GFC_1step)
legend('Actual Data','1-Step Ahead Holt-Winter Forecast with Structural Break', 'FontSize',25, 'Location','northwest')
hold off

title('Plot for Holt-Winter 1-step Forecast with Structural Break for GFC', 'FontSize',30)
ylabel('Volume (m^{3})', 'FontSize',25)


%%
function MSFE_optimised = holtoptimiser_2_6_GFC_1step(v)

Concrete = readmatrix('Concrete_data.csv');
y = Concrete(9:118, 2);
T = length(y);
h1 = 1;
T0 = 8;

syhat_2_6_GFC_1step = zeros(T-h1-T0+1,1);
ytph_2_6_GFC_1step = y(T0+h1:end); % observed y_{t+h}
alpha = v(1); beta = v(2); gamma = v(3); % smoothing parameters
s = 4; % periodicity of seasonality
Lt = mean(y(1:s)); bt = 0; St(1:s) = y(1:s)-Lt;
GFC = 56; % start of GFC period 2008Q4
psi = .7; % parameter that exponentially becomes smaller as period after GFC increases
for t = s+1:T-h1
    if (t < GFC)
	    newLt = alpha*(y(t)-St(t-s)) + (1-alpha)*(Lt+bt); % level
	    newbt = beta*(newLt-Lt) + (1-beta)*bt; % slope
        St(t) = gamma*(y(t)-newLt) + (1-gamma)*St(t-s); % Seasonality
        yhat = newLt + h1*newbt + St(t+h1-s);
	    Lt = newLt; bt = newbt; % update Lt and bt
    elseif (t >= GFC)
        newLt = alpha*(y(t)-St(t-s)) + (1-alpha)*(Lt+bt); % level
	    newbt = beta*(newLt-Lt) + (1-beta)*bt; % slope
        St(t) = gamma*(y(t)-newLt) + (1-gamma)*St(t-s); % Seasonality
        yhat = newLt + h1*newbt + St(t+h1-s) - (15595)*(psi^(t+1-GFC));
	    Lt = newLt; bt = newbt; % update Lt and bt
    end
    
    if t>= T0 % store the forecasts for t >= T0
		syhat_2_6_GFC_1step(t-T0+1,:) = yhat;
 	end
end

MSFE_optimised = mean((ytph_2_6_GFC_1step-syhat_2_6_GFC_1step).^2);

end
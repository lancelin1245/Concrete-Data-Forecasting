%% Forecast for 2023Q1
h3=3;
syhat_2_6_GFC_3step = zeros(T-h3-T0+1,1);
ytph_2_6_GFC_3step = y(T0+h3:end); % observed y_{t+h}
alpha = 0.4391; beta = 0.1636; gamma = 0.2123; % smoothing parameters
s = 4; % periodicity of seasonality
Lt = mean(y(1:s)); bt = 0; St(1:s) = y(1:s)-Lt;
GFC = 56; % start of GFC period 2008Q4
psi = .7; % parameter that exponentially becomes smaller as period after GFC increases
for t = s+1:T
    if (t < GFC)
	    newLt = alpha*(y(t)-St(t-s)) + (1-alpha)*(Lt+bt); % level
	    newbt = beta*(newLt-Lt) + (1-beta)*bt; % slope
        St(t) = gamma*(y(t)-newLt) + (1-gamma)*St(t-s); % Seasonality
        yhat = newLt + h3*newbt + St(t+h3-s);
	    Lt = newLt; bt = newbt; % update Lt and bt
    elseif (t >= GFC)
        newLt = alpha*(y(t)-St(t-s)) + (1-alpha)*(Lt+bt); % level
	    newbt = beta*(newLt-Lt) + (1-beta)*bt; % slope
        St(t) = gamma*(y(t)-newLt) + (1-gamma)*St(t-s); % Seasonality
        yhat = newLt + h3*newbt + St(t+h3-s) - (15595)*(psi^(t+1-GFC));
	    Lt = newLt; bt = newbt; % update Lt and bt
    end
    
    if t>= T0 % store the forecasts for t >= T0
		syhat_2_6_GFC_3step(t-T0+1,:) = yhat;
 	end
end
syhat_2_6_GFC_3step(end)
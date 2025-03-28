clc; close all; clear;

% Generate input signal r(t) (Multi-sine signal)
N = 3000;
M = 1.8;
x = 0.69969;

Band = [0 x];
Range = [-M,M];

SineData = [100,100,1];

% define input signal as multi sine
r = idinput(N, 'prbs', Band, Range);

% Get system response using the pcode file
[u, y] = assignment_sys_33(r, "open loop");

% Create IDDATA object for training data
data = iddata(y, u);

% Set BJ options using bjOptions
options = bjOptions( ...
    'InitialCondition', 'zero', ...
    'Display', 'on', ...
    'Focus', 'prediction');

% Identify the BJ model using the train data and the options
model_BJ = bj(data, [4 3 4 4 1], options);


% Validate Test 1
% define input signal as multi sine
r_val = idinput(N, 'prbs', Band, Range);
[u_val, y_val] = assignment_sys_33(r_val, "open loop");
val_data = iddata(y_val, u_val);
compare(val_data, model_BJ); 



















% % Compute and plot residuals (validation error)
% residuals = y_val - y_pred;  % Compute residuals (error between predicted and actual output)
% 
% % Plot actual vs predicted output
% figure;
% subplot(2, 1, 1);
% plot(y_val, 'b', 'DisplayName', 'Actual Output');
% hold on;
% plot(y_pred.y, 'r', 'DisplayName', 'Predicted Output');
% legend;
% title('BJ Model Validation: Actual vs Predicted');
% 
% % Plot residuals
% subplot(2, 1, 2);
% plot(residuals);
% title('Residuals (Validation Set)');

% % Compute performance metrics
% mse = mean(residuals.^2);  % Mean Squared Error
% disp(['Mean Squared Error (Validation Set): ', num2str(mse)]);
% 
% % Compute R² (coefficient of determination)
% SS_tot = sum((y_val - mean(y_val)).^2);
% SS_res = sum(residuals.^2);
% R2 = 1 - (SS_res / SS_tot);
% disp(['R^2 (Validation Set): ', num2str(R2)]);
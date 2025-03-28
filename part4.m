clc; close all; clear;

% Generate input signal r(t) (PRBS)
N = 3000;
M = 1.8;
x = 0.69969;

Band = [0 x];
Range = [-M, M];

% Define input signal as PRBS
r_prbs = idinput(N, 'prbs', Band, Range);
[u_prbs, y_prbs] = assignment_sys_33(r_prbs, "open loop");
data_prbs = iddata(y_prbs, u_prbs);

% % Define BJ Model
% model_BJ = bj(data_prbs, [4 3 4 3 1]);
% 
% % Define ARMAX
% model_ARMAX = armax(data_prbs, [4 3 3 1]);
% 
% % Define ARX
% model_ARX = arx(data_prbs, [4 3 1]);

% Define OE Changed the orders to 5 4 1 becasue they give the lowest value
model_OE = oe(data_prbs, [5 4 1]);

% % Define ranges for na, nb, nk
% na_range = 1:5;  % Example range for na (AR order)
% nb_range = 1:5;  % Example range for nb (MA order)
% nk_range = 0:2;  % Example range for nk (delay)
% 
% % Initialize variables to store AIC values and model orders
% aic_values = NaN(length(na_range)*length(nb_range)*length(nk_range), 1);
% model_orders = NaN(length(na_range)*length(nb_range)*length(nk_range), 3);  % Store orders [na, nb, nk]
% 
% % Counter for storing the index of the best model
% counter = 1;

% % Loop over all combinations of na, nb, nk
% for na = na_range
%     for nb = nb_range
%         for nk = nk_range
%             % Estimate the OE model with current orders
%             model_OE = oe(data_prbs, [na nb nk]);  % Estimate OE model
% 
%             % Calculate AIC using MATLAB's built-in aic function
%             aic_value = aic(model_OE);  % Compute AIC for the model
% 
%             % Store AIC value and corresponding model orders
%             aic_values(counter) = aic_value;
%             model_orders(counter, :) = [na nb nk];
% 
%             % Increment counter
%             counter = counter + 1;
%         end
%     end
% end

% Find the model with the minimum AIC
% [~, best_idx] = min(aic_values);
% best_orders = model_orders(best_idx, :);
% 
% % Display the best model orders and AIC value
% disp(['Best Model Orders (na, nb, nk): ', num2str(best_orders)]);
% disp(['Best AIC: ', num2str(aic_values(best_idx))]);

% Generate validation input signal
r_val = idinput(N, 'rgs', Band, Range);
[u_val, y_val] = assignment_sys_33(r_val, "open loop");
val_data = iddata(y_val, u_val);

% Validation Test 1 : Compare models with validation data
compare(val_data, model_OE);


% Validation Test 2: Residual Test
resid(val_data, model_OE);


% Minimum Variance Estimate

getcov(model_OE);
% TODO: add reference data.mat and assignment_sys_33.p
% They all should be in a zip file called group_33_ssml.zip


% -------- FINAL REPORT - GROUP 33: SADRA MOOSAVI LAR (2139901) --------- %

% Repository: https://github.com/sadra-hub/5smb0-systemIdentification

% ------------------------------ PART 1 --------------------------------- % 

clc; close all; clear;

% Define an input signal to excite the system
r = -300:0.01:300;
[u, ~] = assignment_sys_33(r,'open loop');
M = max(u);

% Butterworth filter:
numerator = [0.505 1.01 0.505];
denominator = [1 0.7478 0.2722];
F_q = filt(numerator, denominator);

% The bandwidth is the first frequency where 
% the gain drops below 70.79% (-3 dB) of its DC value. (From MATLAB Docs) 
x = bandwidth(F_q)/pi;

figure("Name","PART 1: Bode Plot")
freqz(numerator, denominator)

% Display
disp("M is : " + M)
disp("x is :" + x)

% ------------------------------ PART 2 --------------------------------- %

clc; close all; clear;

% defining datapoints, butterfilter range and bandwidth
N = 1500;
M = 1.8;
x = 0.69969;

Band = [0 x];
Range = [-M,M];

% (at least) 100 frequencies within the bandwidth
SineData = [100,100,1];

% define input signal as multi sine
r = idinput(N,'sine',Band, Range, SineData);

% excite the system using the input signal
[u , y] = assignment_sys_33(r,"open loop");

% create a data object using the input and output of the system
data = iddata(y, u);

% estimate the Frequency Response Function (FRF)
G_frf = spa(data); 

% display the Bode plot of the identified system
figure("Name", "Bode Plot of Identified FRF");
bode(G_frf);
grid on;
title('Bode Plot of the Identified System');

[Sy, f] = cpsd(y, y, [], [], N);    % Output PSD
[Syu, ~] = cpsd(y, u, [], [], N);   % Cross PSD
[Su, ~] = cpsd(u, u, [], [], N);    % Input PSD

% Compute noise spectrum estimate
Phi_v = Sy - abs(Syu).^2 ./ Su;

% Plot the magnitude in a Bode diagram
figure;
semilogx(f, 10*log10(abs(Phi_v)));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Estimated Noise Spectrum \Phi_v');
legend('Noise Spectrum');
grid on;


% ------------------------------ PART 3 --------------------------------- %

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

% Save the iddata object to a .mat file
save('data_open_loop.mat', 'data');

% ------------------------------ PART 4 --------------------------------- %

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

% Define BJ Model
model_BJ = bj(data_prbs, [4 3 4 3 1]);

% Define ARMAX
model_ARMAX = armax(data_prbs, [4 3 3 1]);

% Define ARX
model_ARX = arx(data_prbs, [4 3 1]);

% Define OE 
% Changed the orders to 5 4 1 for open loop
% Changed the orders to 5 5 0 for closed loop
model_OE = oe(data_prbs, [5 5 0]);

% Define ranges for na, nb, nk
na_range = 1:5;  % Example range for na (AR order)
nb_range = 1:5;  % Example range for nb (MA order)
nk_range = 0:2;  % Example range for nk (delay)

% Initialize variables to store AIC values and model orders
aic_values = NaN(length(na_range)*length(nb_range)*length(nk_range), 1);
model_orders = NaN(length(na_range)*length(nb_range)*length(nk_range), 3);  % Store orders [na, nb, nk]

% Counter for storing the index of the best model
counter = 1;

% Loop over all combinations of na, nb, nk
for na = na_range
    for nb = nb_range
        for nk = nk_range
            % Estimate the OE model with current orders
            model_OE = oe(data_prbs, [na nb nk]);  % Estimate OE model

            % Calculate AIC using MATLAB's built-in aic function
            aic_value = aic(model_OE);  % Compute AIC for the model

            % Store AIC value and corresponding model orders
            aic_values(counter) = aic_value;
            model_orders(counter, :) = [na nb nk];

            % Increment counter
            counter = counter + 1;
        end
    end
end

% Find the model with the minimum AIC
[~, best_idx] = min(aic_values);
best_orders = model_orders(best_idx, :);

% Display the best model orders and AIC value
disp(['Best Model Orders (na, nb, nk): ', num2str(best_orders)]);
disp(['Best AIC: ', num2str(aic_values(best_idx))]);

% Generate validation input signal
r_val = idinput(N, 'rgs', Band, Range);
[u_val, y_val] = assignment_sys_33(r_val, "open loop");
val_data = iddata(y_val, u_val);

% Validation Test 1 : Compare models with validation data
compare(val_data, model_OE);


% Validation Test 2: Residual Test
% resid(val_data, model_OE);


% Minimum Variance Estimate

getcov(model_OE);


% ------------------------------ PART 5 --------------------------------- %

clc; close all; clear;

% Initialize variables for Monte Carlo Simulation
N = 3000;               % Number of data points
M = 1.8;                % Maximum value of input signal
x = 0.69969;            % Bandwidth range

Band = [0 x];          % Frequency range
Range = [-M, M];       % Range of values for input signal
num_simulations = 1;   % Number of Monte Carlo simulations

% Define the model orders for the OE model
nb = 5;  % Order of the numerator polynomial
nf = 4;  % Order of the denominator polynomial
nk = 1;  % Input-output delay

% Initialize arrays to store parameter estimates from each run
B_estimates = zeros(num_simulations, nb);
F_estimates = zeros(num_simulations, nf);

% Perform Monte Carlo simulations
for i = 1:num_simulations

    % Define input signal as PRBS
    r_prbs = idinput(N, 'prbs', Band, Range);
    [u_prbs, y_prbs] = assignment_sys_33(r_prbs, "open loop");
    data_prbs = iddata(y_prbs, u_prbs);
    
    % Estimate the OE model for the data
    model_OE = oe(data_prbs, [nb, nf, nk]);
    
    % Store the estimated parameters
    B_estimates(i, :) = model_OE.B(2:end);  % Numerator coefficients
    F_estimates(i, :) = model_OE.F(2:end);  % Denominator coefficients
end

% Calculate the mean and variance of the estimated parameters
B_mean = mean(B_estimates, 1);
B_variance = var(B_estimates, 0, 1);
F_mean = mean(F_estimates, 1);
F_variance = var(F_estimates, 0, 1);

% Covariance of matrix
cov_matrix = getcov(model_OE);

% Display results
disp('Mean of B coefficients:');
disp(B_mean);
disp('Variance of B coefficients:');
disp(B_variance);

disp('Mean of F coefficients:');
disp(F_mean);
disp('Variance of F coefficients:');
disp(F_variance);

disp('Covarince of Model OE:')
disp(cov_matrix)

% The diagonal of the covariance matrix gives the parameter variances
parameter_variances = diag(cov_matrix);
disp('Estimated Parameter Variances:');
disp(parameter_variances);

% ------------------------------ PART 6 --------------------------------- %

clc; close all; clear;

% Define the parameters
N = 3000;               % Number of data points
M = 1.8;                % Maximum value of input signal
x = 0.69969;            % Bandwidth range

Band = [0 x];           % Frequency range
Range = [-M, M];        % Range of values for input signal

% Generate a random reference signal for the closed-loop system
r_prbs = idinput(N, 'prbs', Band, Range);

% Get closed-loop and open-loop data using the assignment_sys_33 function
[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'closed loop');
data_prbs_close = iddata(y_prbs, u_prbs); 

[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'open loop');
data_prbs_open = iddata(y_prbs, u_prbs);

% perform system identification for the closed-loop system using the OE model
% model orders are optimized by AIC value

model_OE_close = oe(data_prbs_close, [5 5 0]);
model_OE_open = oe(data_prbs_open, [5 4 1]);

% Compare the closed-loop and open-loop models against validation data
r_prbs = idinput(N, 'prbs', Band, Range);
[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'closed loop');
data_prbs__close_val = iddata(y_prbs, u_prbs); 

[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'open loop');
data_prbs__open_val = iddata(y_prbs, u_prbs); 


subplot(2, 1, 1);
compare(data_prbs__close_val, model_OE_open, model_OE_close);
legend('Closed-Loop Data', 'Open-Loop Model', 'Closed-Loop Model');

subplot(2, 1, 2);
compare(data_prbs__open_val, model_OE_open, model_OE_close);
legend('Open-Loop Data', 'Open-Loop Model', 'Closed-Loop Model');

% Display results of the closed-loop identification
disp('Closed-Loop Model Parameters:');
disp(model_OE_close);
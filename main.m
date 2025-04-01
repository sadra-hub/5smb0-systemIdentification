% TODO: add reference data.mat and assignment_sys_33.p
% They all should be in a zip file called group_33_ssml.zip


%% -------- FINAL REPORT - GROUP 33: SADRA MOOSAVI LAR (2139901) -------- %

% Repository: https://github.com/sadra-hub/5smb0-systemIdentification

%% ------------------------------ PART 1 -------------------------------- % 

clc; close all; clear;

% Define an input signal to excite the system
r = linspace(-100,100,1000);
[u, ~] = assignment_sys_33(r,'open loop');

% Butterworth filter:
numerator = [0.505 1.01 0.505];
denominator = [1 0.7478 0.2722];
F_q = filt(numerator, denominator);

% M is the maximum or minmum response 
M = max(u);

% The bandwidth is the first frequency where 
% the gain drops below 70.79% (-3 dB) of its DC value. 
% _From MATLAB Docs
x = bandwidth(F_q)/pi;

figure("Name", "Bode Plot of F(q)");
title('$F(q)$', 'Interpreter', 'latex');
freqz(numerator, denominator);

% Display
disp("M is : " + M)
disp("x is :" + x)

%% ------------------------------ PART 2 -------------------------------- %

% clc; close all; clear;

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


%% ------------------------------ PART 3 -------------------------------- %

% clc; close all; clear;

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
% data is saved to a NEW file, to not overwrite the old data which is used 
% for figures and report
save('data/data_open_loop_NEW.mat', 'data');

%% ------------------------------ PART 4 -------------------------------- %

% clc; close all; clear;

% Generate input signal r(t) (PRBS)
N = 3000;
M = 1.8;
x = 0.69969;

Band = [0 x];
Range = [-M, M];

% Define input signal as PRBS (loading the produced data in part 2)
data_open_loop = load('data/data_open_loop.mat');
data_prbs = data_open_loop.data;

% % Define search ranges for BJ model orders
% nb_range = 1:4;
% nc_range = 1:5;
% nd_range = 1:5;
% nf_range = 1:4;
% nk_range = 0:2;
% 
% % Preallocate storage
% total_combinations = length(nb_range) * length(nc_range) * length(nd_range) * length(nf_range) * length(nk_range);
% aic_values = NaN(total_combinations, 1);
% model_orders = NaN(total_combinations, 5);  % [nb nc nd nf nk]
% 
% % Counter
% counter = 1;
% 
% % Grid search over all combinations
% for nb = nb_range
%     for nc = nc_range
%         for nd = nd_range
%             for nf = nf_range
%                 for nk = nk_range
%                     try
%                         model_BJ = bj(data_prbs, [nb nc nd nf nk]);
%                         aic_val = aic(model_BJ);
% 
%                         % Store values
%                         aic_values(counter) = aic_val;
%                         model_orders(counter, :) = [nb nc nd nf nk];
%                     catch
%                         % If model fails, skip and continue
%                         aic_values(counter) = inf;
%                     end
%                     counter = counter + 1;
%                 end
%             end
%         end
%     end
% end

% best orders according to the code above which I commented out 
% to save time
best_orders = [4,5,5,4,1];

% Display result
disp('Best BJ Model Orders [nb nc nd nf nk]:');
disp(best_orders);

% fit final best model
model_BJ = bj(data_prbs, best_orders);

% load open loop validation data for validitation tests
% validation data are PRBS signal similar to input
data_open_loop_validation = load('data/data_open_loop_validation.mat');
val_data_open = data_open_loop_validation.val_data;

% Validation Test 1 : Pole 
figure;
pzmap(model_BJ);
title('Pole-Zero Map of BJ Model');
grid on;

is_stable = isstable(model_BJ);
disp(['Is the BJ model stable? ', mat2str(is_stable)]);

% Validation Test 2: Residual Analysis 
figure;
resid(val_data_open, model_BJ);
sgtitle('Residual Analysis of BJ Model');

% Compare with results in Part 2
figure;
bode(G_frf, model_BJ);
grid on;
legend('Nonparametric FRF (spa)', 'Parametric BJ Model');
title('Comparison of Nonparametric FRF and BJ Model');

% Extract parameter variances (diagonal elements)
cov_matrix = getcov(model_BJ);
param_variances = diag(cov_matrix);
disp('Parameter variances (theoretical):');
disp(param_variances);


%% ------------------------------ PART 5 -------------------------------- %

% clc; close all; clear;

% Initialize variables for Monte Carlo Simulation
N = 3000;               % Number of data points
M = 1.8;                % Saturation limit
x = 0.69969;            % Butterworth bandwidth (normalized)

Band = [0 x];           
Range = [-M, M];        
num_simulations = 100;   % Number of Monte Carlo runs

% Define best BJ model orders (from previous AIC search)
nb = 4;  % B(q) numerator order
nc = 5;  % C(q) noise numerator order
nd = 5;  % D(q) noise denominator order
nf = 4;  % F(q) system denominator order
nk = 1;  % Input-output delay

% Get total number of parameters
num_BJ_params = nb + nc + nd + nf;

% Initialize array to store parameter estimates
BJ_parameters = zeros(num_simulations, num_BJ_params);

% Monte Carlo simulation
for i = 1:num_simulations
    r_prbs = idinput(N, 'prbs', Band, Range);
    [u, y] = assignment_sys_33(r_prbs, 'open loop');
    data = iddata(y, u);
    
    try
        model_BJ = bj(data, [nb nc nd nf nk]);
        BJ_parameters(i, :) = model_BJ.ParameterVector';  % Store all parameters
    catch
        warning(['Model estimation failed at iteration ', num2str(i)]);
        BJ_parameters(i, :) = NaN;  % Handle occasional fitting failures
    end
end

% Remove failed runs (rows with NaNs)
BJ_parameters = BJ_parameters(~any(isnan(BJ_parameters), 2), :);

% Compute mean and variance of estimated parameters
BJ_mean = mean(BJ_parameters, 1);
BJ_variance = var(BJ_parameters, 0, 1);

% Theoretical covariance matrix (from one final fit)
model_BJ = bj(data, [nb nc nd nf nk]);
cov_matrix = getcov(model_BJ);
parameter_variances = diag(cov_matrix)';

% Extract indices for B and F parameters (per question 5.3)
% B starts at index 1, F starts after B + C + D
B_indices = 1:nb;
F_indices = (nb + nc + nd + 1):(nb + nc + nd + nf);

% Extract Monte Carlo and theoretical variances for B and F
B_variance_mc = BJ_variance(B_indices);
F_variance_mc = BJ_variance(F_indices);

B_variance_theoretical = parameter_variances(B_indices);
F_variance_theoretical = parameter_variances(F_indices);

% Combine for plotting
all_variances_mc = [B_variance_mc, F_variance_mc];
all_variances_theoretical = [B_variance_theoretical, F_variance_theoretical];

% Create bar plot
figure;
bar(1:(nb+nf), [all_variances_mc(:), all_variances_theoretical(:)]);
xlabel('Parameter Index for B(q) & F(q)');
ylabel('Variance');
title('Monte Carlo vs. Theoretical Variance');
legend('Monte Carlo Variance', 'Theoretical Variance');
grid on;

%% ------------------------------ PART 6 -------------------------------- %

% clc; close all; clear;

% load closed loop data (this is a PRBS signal)
data_closed_loop = load('data/data_closed_loop.mat');
data_prbs_close = data_closed_loop.data;

% perform system identification for the closed-loop system
% model orders are optimized by minimum AIC value
model_ARMAX = armax(data_prbs_close, [5 5 30 1]);

disp(aic(model_ARMAX))

% load validation data for closed loop
data_closed_loop_validation = load('data/data_closed_loop_validation.mat');
val_data_closed = data_closed_loop_validation.val_data;


% Validation Test 1 : Pole 
figure;
pzmap(model_ARMAX);
title('Pole-Zero Map of ARMAX Model');
grid on;

is_stable = isstable(model_ARMAX);
disp(['Is the ARMAX model stable? ', mat2str(is_stable)]);

% Validation Test 2: Residual Analysis 
figure;
resid(val_data_closed, model_ARMAX);
sgtitle('Residual Analysis of ARMAX Model');

% estimate the Frequency Response Function (FRF)
G_frf_closed = spa(data_prbs_close);

figure;
bode(G_frf_closed, model_ARMAX);
legend('Nonparametric FRF (spa)', 'Parametric ARMAX Model');
title('Comparison of Nonparametric FRF and ARMAX Model');
grid on;

figure;
compare(val_data_closed, model_ARMAX);
title('ARMAX Model Performance on New Validation Data');
legend('Measured Output', 'Model Prediction');

figure;
compare(val_data_open, model_BJ);
title('BJ Model Performance on New Validation Data');
legend('Measured Output', 'Model Prediction');
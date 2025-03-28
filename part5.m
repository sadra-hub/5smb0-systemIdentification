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

clc; close all; clear;

% Define the parameters
N = 3000;               % Number of data points
M = 1.8;                % Maximum value of input signal
x = 0.69969;            % Bandwidth range

Band = [0 x];          % Frequency range
Range = [-M, M];       % Range of values for input signal

% Generate a random reference signal for the closed-loop system
r_prbs = idinput(N, 'prbs', Band, Range);  % Use PRBS for reference signal

% Get closed-loop data using the modified function
[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'closed loop');  % 'closed loop' mode
data_prbs_close = iddata(y_prbs, u_prbs);  % Create the closed-loop data object

[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'open loop');
data_prbs_open = iddata(y_prbs, u_prbs);

% perform system identification for the closed-loop system using the OE model
model_OE_close = oe(data_prbs_close, [5 5 0]);
model_OE_open = oe(data_prbs_open, [5 4 1]);

% Compare the closed-loop and open-loop models against new data
r_prbs = idinput(N, 'prbs', Band, Range);
[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'closed loop');  % 'closed loop' mode
data_prbs__close_val = iddata(y_prbs, u_prbs);  % Create the closed-loop data object

[u_prbs, y_prbs] = assignment_sys_33(r_prbs, 'open loop');  % 'closed loop' mode
data_prbs__open_val = iddata(y_prbs, u_prbs);  % Create the closed-loop data object


subplot(2, 1, 1);
compare(data_prbs__close_val, model_OE_open, model_OE_close);
legend('Closed-Loop Data', 'Open-Loop Model', 'Closed-Loop Model');

subplot(2, 1, 2);
compare(data_prbs__open_val, model_OE_open, model_OE_close);
legend('Open-Loop Data', 'Open-Loop Model', 'Closed-Loop Model');

% Display results of the closed-loop identification
disp('Closed-Loop Model Parameters:');
disp(model_OE_close);
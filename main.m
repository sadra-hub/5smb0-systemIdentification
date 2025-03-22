% TODO: add reference data.mat and assignment_sys_33.p
% They all should be in a zip file called group_33_ssml.zip

close all
clear
clc

% 1.1 DEFINING F(q)
% Define coefficients of the transfer function
num = [0.505, 1.01, 0.505];  % Numerator coefficients (for q^-2, q^-1, and q^0)
den = [1, 0.7478, 0.2722];   % Denominator coefficients (for q^-2, q^-1, and q^0)

% Define the transfer function F(q)
F = tf(num, den, 1); % '1' indicates it's a discrete transfer function

% Plot Bode diagram
% figure;
% bode(F);
% grid on;
% title('Bode Diagram of F(q)');

% Find the cutoff frequency (where the gain is -3 dB)
[mag, phase, wout] = bode(F);
mag = squeeze(mag); % Magnitude at different frequencies

% Convert magnitude to dB and find the index where the magnitude is -3 dB
mag_dB = 20 * log10(mag); 
cutoff_idx = find(mag_dB <= -3, 1, 'first'); % First index where magnitude drops to -3 dB

% Display the cutoff frequency (in terms of normalized frequency)
cutoff_freq = wout(cutoff_idx) / pi;  % Normalized frequency in terms of pi

fprintf('The cutoff frequency is %.3fπ.\n', cutoff_freq);


% 1.2 FINDING M

% Define the simulation parameters
T_sim = 50;  % Total simulation time (seconds)
time = 0:T_sim;  % Time vector
M_initial_guess = 2;  % Initial guess for the saturation limit (M)

% Define the ramp signal r(t) for testing
r = @(t) 0.2 * t;  % Ramp signal with a slope of 0.2

% Pre-allocate output vectors
u = zeros(size(time));  % System output after saturation : u(t)

% Simulate the system and apply saturation
for i = 1:length(time)
    % Call the system function with r(t) at each time step
    [u_temp, ~] = assignment_sys_33(r(time(i)), 'open loop');
    
    % Apply the saturation function S(x) to the output
    u(i) = min(max(u_temp, -M_initial_guess), M_initial_guess);
end

% Plot the output to observe the saturation effect
figure;
plot(time, u);
xlabel('Time (s)');
ylabel('Actuator Output (u(t))');
title('System Response with Saturation');

% Find the time at which the output saturates
[max_output, max_index] = max(u);
fprintf('Saturation occurs around time %.2f seconds with output %.2f\n', time(max_index), max_output);

% Estimate the saturation limit M based on the output behavior
% Assuming the output value when it stops increasing corresponds to M
M_estimated = max_output;
fprintf('Estimated saturation limit M = %.2f\n', M_estimated);


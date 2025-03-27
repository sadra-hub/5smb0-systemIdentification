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
% figure("Name", "Bode Plot of Identified FRF");
% bode(G_frf);
% grid on;
% title('Bode Plot of the Identified System');

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
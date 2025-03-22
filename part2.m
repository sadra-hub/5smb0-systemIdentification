clc 
close all
clear

% Generate multi-sine input signal
N = 1500;
t = (0:N-1)'; % Time vector
% the passband here is below the cut off frequency we calculated at part 1
frequencies = linspace(0.01, 0.719, 100); % 100 frequencies in passband
r = sum(sin(2 * pi * frequencies .* t), 2); % Multi-sine signal

% Obtain system response
[u, y] = assignment_sys_33(r, "open loop");

% Estimate FRF using nonparametric identification
[H, f] = tfestimate(u, y, [], [], [], 1);

% % Bode plot
% figure;
% subplot(2,1,1);
% semilogx(f, 20*log10(abs(H))); % Magnitude response
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('Bode Plot - Magnitude Response');
% grid on;
% 
% subplot(2,1,2);
% semilogx(f, angle(H)*180/pi); % Phase response
% xlabel('Frequency (Hz)');
% ylabel('Phase (degrees)');
% title('Bode Plot - Phase Response');
% grid on;

[Pyy, f] = cpsd(y, y, [], [], N);
[Pyu, ~] = cpsd(y, u, [], [], N);
[Puu, ~] = cpsd(u, u, [], [], N);

% Compute the noise spectrum:
Pvv = Pyy - abs(Pyu).^2 ./ Puu;

% TODO : Why it mentions the welch avergaing mehtod? am I missing somethign
% here? 
% visualize the noise spectrum
figure;
loglog(f, abs(Pvv));
xlabel('Frequency (Hz)');
ylabel('Magnitude of \Phi_v(\omega)');
title('Estimated Noise Spectrum');
grid on;


clc; close all; clear;

% Defina an input signal to excite the system
r = -300:0.01:300;
[u, y] = assignment_sys_33(r,'open loop');
M = max(u);

% Butterworth filter:
numerator = [0.505 1.01 0.505];
denominator = [1 0.7478 0.2722];
F_q = filt(numerator, denominator);

% The bandwidth is the first frequency where 
% the gain drops below 70.79% (-3 dB) of its DC value. (From MATLAB) 
x = bandwidth(F_q)/pi;

figure("Name","Bode Plot")
freqz(numerator, denominator)

% Display
disp("M is : " + M)
disp("x is :" + x)

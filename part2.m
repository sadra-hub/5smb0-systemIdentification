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




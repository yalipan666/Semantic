function [y] = sub_notch(x,fs,f_notch,Q_factor)
% This function is used to design a second-order IIR notch filter 
% x: input
% fs: sampling rate
% f_notch: the frequency to be removed
% Q_factor: to adjust the bandwidth 'BW' of the filter
% y: filtered signal

% Design a notch filter with a bandwidth of BW at a level -Ab in decibels to remove a f_notch-Hz tone
Wo = f_notch/(fs/2)*pi; % the normalized frequency 
BW = Wo/Q_factor;       % the bandwidth of the filter
Ab = -3;      % -3dB width
Gb   = 10^(Ab/20);      % filter design
beta = (sqrt(1-Gb.^2)/Gb)*tan(BW/2);
gain = 1/(1+beta);
b  = gain*[1 -2*cos(Wo) 1];             % the numerator coefficient vector b
a  = [1 -2*gain*cos(Wo) (2*gain-1)];    % the denominator coefficient vector a
nb = length(b)-1;       % the length of b
na = length(a)-1;       % the length of a

% IIR filtering
y = filtfilt(b,a,x);


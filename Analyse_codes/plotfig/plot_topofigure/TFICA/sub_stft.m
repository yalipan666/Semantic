function [S, U, P, F] = sub_stft(x, xtimes, t, f, Fs, winsize, detrend_opt, wintype, padvalue)
% This sub-function is used to estimate time-varying complex spectrum using 
% short-time Fourier transform (STFT)
% Note that it is only used for multi-channel EEG data.

% /Input/
% x: the original data samples (Time Points x Channels)
% xtimes: the time axis of the original data
% t: evaluated time points in STFT
% f: evaluated frequency bins in STFT
% Fs: sampling rate
% winsize: window size (NOTE: the unit is sec)
% detrend_opt: detrent or not ('1' for detrend is default, '0' for not)
% wintype: window type (default: hanning window)
% padvalue: the pad value used for padding data (default: '0'; See matlab function 'padarray')

% /Output/
% S: complex time-frequency value
% U: a constant to normalize the power
% P: squared magnitude without phase
% F: phase information


% fprintf('\nShort-time Fourier Transform: ')

%% Parameters
if nargin<=6; detrend_opt = 1; end      % detrend (1, default) the data or not (0)
if nargin<=7; wintype = 'hann'; end;    % default window type is 'hann'
if nargin<=8; padvalue = 0; end;        % default padding value is 0

if size(x,1)==1; x = x.'; end; % transpose data if the 1st dimension (time) is 1
N_Trials = size(x,2); % number of trials
L = round(winsize*Fs/2); % half of window size (points)
h = L*2+1; %  window size (points)
win = window(wintype,h); % window (one trial)
if ismember(0,win); win = window(wintype,h+2); win = win(find(win~=0)); end % remove zeros points in the window
W = repmat(win,1,N_Trials); % window (all trials)
U = win'*win;  % compensates for the power of the window

% index of time points in time-frequency domain
dt = t(2)-t(1); % time interval (uniform step)
[junk,t_idx_min] = min(abs(xtimes-t(1)));
[junk,t_idx_max] = min(abs(xtimes-t(end)));
t_idx = t_idx_min:round(dt*Fs):t_idx_max; 
N_T = length(t);

% index of time points in time-frequency domain
df = f(2)-f(1); % frequency step (uniform step)
nfft = round(Fs/df) * max(1,2^(nextpow2(h/round(Fs/df)))); % points of FFT
f_full = Fs/2*linspace(0,1,round(nfft/2)+1);
[junk,f_idx_min] = min(abs(f_full-f(1)));
[junk,f_idx_max] = min(abs(f_full-f(end)));
f_idx = f_idx_min:max(1,2^(nextpow2(h/round(Fs/df)))):f_idx_max; 
if numel(f_idx)>numel(f); f_idx = f_idx(1:end-1); end
N_F = length(f);

% fprintf('%d Time Points x %d Frequency Bins x %d Trials\n',N_T,N_F,N_Trials);

%% Pad x (default mode is "zero")
if detrend_opt; x = detrend(x); end; % Remove linear trends
X = padarray(x, L, padvalue);        % padding data   
if detrend_opt; X = detrend(X); end; % Remove linear trends

%% STFT
S = single(zeros(N_F,N_T,N_Trials));
% fprintf('Processing...     ')
for n=1:N_T
%     fprintf('\b\b\b\b%3.0f%%',n/N_T*100)
    X_n = X(t_idx(n)+[0:h-1],:); 
	if detrend_opt; X_n = detrend(X_n); end; % Remove linear trends
    S_n = fft(X_n.*W,nfft,1);
    S(:,n,:) = S_n(f_idx,:);
end
S2 = S.*conj(S);
P = S2/Fs/U;
F = S./sqrt(S2);
% fprintf('  Done!\n')
%     P(2:end-1,:,:) = 2*P(2:end-1,:,:); % *2 means to consider the components
%     in negative freq bands. If it is required to let the program is the
%     same as spectrogram of Matlab, then enable it. Note that the end
%     frequency should be the Nyquist frequency


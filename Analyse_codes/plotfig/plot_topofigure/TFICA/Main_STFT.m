clear all
%% Parameters 
SubNames = {'AN';'BS';'DP';'DT';'GDI';'LH';'LT';'SB';'TC';'VD'};
StrConds = {'PAIN','TOUCH'};
N_Subs = numel(SubNames);
N_Conds = numel(StrConds);

for n_sub=1:N_Subs
for n_cond=1:N_Conds

filename = ['Results\EEG\EEG_',SubNames{n_sub},'_',StrConds{n_cond},'.mat'];
load (filename)

%% Set parameters for time-frequency analysis
f_min = 1; % the minimum frequency evaluated in TFA (unit: Hz; it can be altered)
f_max = 60; % the maximum frequency evaluated in TFA (unit: Hz; it MUST be smaller than Fs/2)
df = 1; % the frequency step between two frequency bins evaluated in TFA
f = [f_min:df:min(Fs/2,f_max)]'; % frequency bins evaluated in TFA (it is better to use a frequency step of 1/(2^k), athough others are OK)
N_Freq = numel(f); % number of freqeuncy bins

detrend_opt = 1;    % detrend the data when performing STFT (if reconstruct is needed, then it is 0)
wintype = 'hann';   % hanning window 
padvalue = 0;       % padding data with 0
winsize = 0.2; % <========= window size (unit: sec) used for short-time Fourier transform

[N_Channels,N_Time,N_Trials] = size(x);

%% Time-frequency analysis (STFT)
S = single(zeros(N_Channels,N_Freq,N_Time,N_Trials));
for n_chan=1:N_Channels
fprintf('\nChannel: %s (%2.0f of %2.0f)\n',ChannelLabels{n_chan},n_chan,N_Channels)
[S(n_chan,:,:,:),U] = sub_stft(squeeze(x(n_chan,:,:)), xtimes, t, f, Fs, winsize, detrend_opt); % S is the complex time-frequency data
end % end of channel

%% Save TFD
savedir = 'Results\STFT\'; % LR for low resolution (2Hz step)
if (exist(savedir)==0); mkdir (savedir); end;
savefile = [savedir,'STFT_',SubNames{(n_sub)},'_',StrConds{n_cond}];
save (savefile,'S','U','xtimes','t','f','Fs','winsize','ChannelLabels');

end % end of condition
end % end of selected subjects
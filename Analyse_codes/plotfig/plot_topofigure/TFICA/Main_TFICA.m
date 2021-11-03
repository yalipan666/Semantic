clear all
% TF-ICA

%% Experimental parameters
SubNames = {'VD'};
StrConds = {'TOUCH'};
N_Subs = numel(SubNames);
N_Conds = numel(StrConds);

DirName = ['Results\'];

TF_Method = 'STFT';

%% LOAD DATA
for n_cond=1:N_Conds 
for n_sub=1:N_Subs
eegfile = [DirName,'EEG\EEG_',SubNames{(n_sub)},'_',StrConds{n_cond}]; % <========= Select a file
load (eegfile)
tfdfile = [DirName,TF_Method,'\',TF_Method,'_',SubNames{(n_sub)},'_',StrConds{n_cond}] % <========= Select a file
load (tfdfile)
% S: chan x freq x time x trial; TIME_FREQUENCY DATA BY STFT
% x: chan x time x trial;        EEG DATA
% t: time axis (sec); xtimes = t;
% f: frequency bin (Hz)
% winsize: Window size used in STFT
% ChannelLabels: channel labels
% Fs: sampling rate

% truncate data (if needed)
if size(t,1)==1; t = t.'; end;
if size(xtimes,1)==1; xtimes = xtimes.'; end;
t_lim = [-1,1];
f_lim = [1,Fs/2];
t_idx = find(t>=t_lim(1) & t<=t_lim(2));
f_idx = find(f>=f_lim(1) & f<=f_lim(2));
S = single(S(:,f_idx,t_idx,:));
x = single(x(:,t_idx,:));
t = t(t_idx); xtimes = xtimes(t_idx);

% parameters
[N_Channels, N_Freq, N_Time, N_Trials] = size(S);
% for TFICA
complexmixing= 0;    % Do we use complex-valued mixing (instead of real-valued)? '1'-yes,'0'-no
removeoutliers= 0;   % Another option: Remove outliers? '1'-yes,'0'-no
maxiter=1000;        % Maximum number of iterations in estimation loop

% for STFT
detrend_opt = 1;    % detrend the data when performing STFT
wintype = 'hann';   % hanning window 
padvalue = 0;       % padding data with 0

%% METHOD 1: Infomax-ICA
% dimension of data used in runica: (chans, frames, epochs)
% remove the baseline
dataout = rmbase(x,size(x,2),[]);
NumICs = N_Channels;
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(dataout,'pca',NumICs);
% re-organize the outputs
s_T = single(reshape(activations,[NumICs,N_Time,N_Trials]));
clear activations dataout
% In runica: activations = weights*sphere*raw_data
A_T = single(pinv(weights*sphere));
W_T = single(weights);
% STFT of ICA-sources
S_T = single(zeros(NumICs,N_Freq,N_Time,N_Trials)); % TFD data
for n_ic=1:NumICs
    S_T(n_ic,:,:,:) = sub_stft(squeeze(s_T(n_ic,:,:)), xtimes, t, f, Fs, winsize, detrend_opt); % S is the complex time-frequency data
end


%% METHOD 2: Fourier-ICA
pcadim=N_Channels;   % Dimension after PCA on channels
NumICs=pcadim;   % Number of components in decomposition
[S_TF,A_TF,W_TF] = sub_fourierica(S, NumICs, pcadim, complexmixing, removeoutliers, maxiter);
% dimension of A is (chan x IC); dimension of W is (IC x chan);
A_TF = single(A_TF); W_TF = single(W_TF);
% inverse STFT
s_TF = single(zeros(NumICs,N_Time,N_Trials));
for n_ic=1:NumICs
    s_TF(n_ic,:,:) = sub_istft(squeeze(S_TF(n_ic,:,:,:)), xtimes, t, f, Fs, winsize);
end


%% METHOD 3: Probabilistic Infomax-ICA
% Selection of Number of InfomaxICs
prob_t=sub_pca_dim(x);
% tmpdata is EEG data (2 dimension: channel num*data length) data length= epoch length*epoch num
[tmp,NumIC_T]=max(prob_t.rrn); % IC_num is the estimated number of ICs
% dimension of data used in runica: (chans, frames, epochs)
% remove the baseline
dataout = rmbase(x,size(x,2),[]);
% ICA
NumICs = NumIC_T;
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(dataout,'pca',NumICs);
% re-organize the outputs
s_PT = single(reshape(activations,[NumICs,N_Time,N_Trials]));
clear activations dataout
% In runica: activations = weights*sphere*raw_data
A_PT = single(pinv(weights*sphere));
W_PT = single(weights);
% STFT of ICA-sources
S_PT = single(zeros(NumICs,N_Freq,N_Time,N_Trials)); % TFD data
for n_ic=1:NumICs
    S_PT(n_ic,:,:,:) = sub_stft(squeeze(s_PT(n_ic,:,:)), xtimes, t, f, Fs, winsize, detrend_opt); % S is the complex time-frequency data
end


%% METHOD 4: Probabilistic Fourier-ICA
% Selection of Number of ICs
prob_tf=sub_pca_dim(S);
% tmpdata is EEG data (2 dimension: channel num*data length) data length= epoch length*epoch num
[tmp,NumIC_TF]=max(prob_tf.rrn); % IC_num is the estimated number of ICs
% NumIC_TF = NumIC_T;
% Probabilistic TF-ICA
pcadim=NumIC_TF;   % Dimension after PCA on channels
NumICs=pcadim;   % Number of components in decomposition
[S_PTF,A_PTF,W_PTF] = sub_fourierica(S, NumICs, pcadim, complexmixing, removeoutliers, maxiter);
% dimension of A is (chan x IC); dimension of W is (IC x chan);
A_PTF = single(A_PTF); W_PTF = single(W_PTF);
% inverse STFT
s_PTF = single(zeros(NumICs,N_Time,N_Trials));
for n_ic=1:NumICs
    s_PTF(n_ic,:,:) = sub_istft(squeeze(S_PTF(n_ic,:,:,:)), xtimes, t, f, Fs, winsize);
end

%% SAVE RESULTS
savedir = ['Results\AllICA\'];
if (exist(savedir)==0); mkdir (savedir); end;
savefile = [savedir,'AllICA_',SubNames{n_sub},'_',StrConds{n_cond},'.mat'];
save (savefile,'S_*','s_*','A_*','W_*','prob*','U','t','f','xtimes','Fs','winsize','ChannelLabels')
clear S_* s_* A_* W_*

%%
end % sub
end % cond

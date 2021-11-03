clear all

%% Parameters 
SubNames = {'AN';'BS';'DP';'DT';'GDI';'LH';'LT';'SB';'TC';'VD'};
StrConds = {'PAIN','TOUCH'};
N_Subs = numel(SubNames);
N_Conds = numel(StrConds);

for n_sub=1:N_Subs
for n_cond=1:N_Conds

filename = ['Data\',SubNames{n_sub},'_',StrConds{n_cond},'2_pre-probe_L.mat']
load (filename)
varname = filename(6:end-4); varname(find(varname=='-')) = '_';
eval(['ERP = ',varname,';']); eval(['clear ',varname,';']); 

%% EEG Parameters
Fs = round(1/ERP.XStep);
xtimes = [1:ERP.XSize]./Fs + ERP.XStart;
x_idx = find(xtimes>=-1 & xtimes<=1);
xtimes = xtimes(x_idx);
N_Channels = numel(ERP.Channels);
ChannelLabels = []; x = [];
for n_chan=1:N_Channels
    ChannelLabels{n_chan,1} = ERP.Channels(n_chan).label;
    x(n_chan,:,:) = squeeze(ERP.RealData(x_idx,1,n_chan,:));
end

%% Decimate data
DS = 8; % 1024/8 = 128
Fs = Fs/DS;
x0 = [];
for n_chan=1:N_Channels
    for n_trial=1:size(x,3)
        eval(['x0(n_chan,:,n_trial) = sub_decimate(squeeze(x(n_chan,:,n_trial)),DS);'])
    end
end
x = x0; clear x0;
xtimes = xtimes(1:DS:end);
t = xtimes;
[N_Channels,N_Time,N_Trials] = size(x);

%% baseline correction
t_baseln_lim = [-0.5, 0]; % <===================================!!!
t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2)));
x = x - repmat(mean(x(:,t_baseln_idx,:),2),[1,size(x,2),1]);

%% high-pass filter
f_cutoff = 0.3;
Wp = f_cutoff/Fs*2;
[b,a]=butter(6,Wp,'high');
for n_chan=1:N_Channels
    for n_trial=1:N_Trials
        x(n_chan,:,n_trial) = single(filtfilt(b,a,double(x(n_chan,:,n_trial))));
    end
end

%% Save EEG
savedir = 'Results\EEG\';
if (exist(savedir)==0); mkdir (savedir); end;
savefile = [savedir,'EEG_',SubNames{n_sub},'_',StrConds{n_cond}];
save (savefile,'x','xtimes','t','Fs','ChannelLabels');

end % end of condition
end % end of selected subjects
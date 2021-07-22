%%% connect dataset that are broken into two parts accidentally
%%% copy from Lexical\ConnectUnexpectedBrokenData

clear

%%% setting parameters
subname = '20210708_b5f6';
megfile = ['of_' subname(10:end)];
savepath = ['Z:\Semantic\Analyse_data\of_' subname '\'];
ptbpath = 'Z:\Semantic\RawData\PTB_data\';
megpath = 'Z:\Semantic\RawData\MEG_data\';

%% %===== connecting the ptb data
%%% load the second part ptb data
load([ptbpath subname '_2.mat'])
broke_trlid = find(~isnan(Result.FixationON),1);
disp(['*** broken trl id: ' num2str(broke_trlid) '***']);
Para_all = Para;
Result_all = Result;

%%% load the first part ptb data
load([ptbpath subname '.mat'])
idx_1 = 1:broke_trlid-1;
Para_all.FixDuration(idx_1) = Para.FixDuration(idx_1);
Result_all.FixationON(idx_1) = Result.FixationON(idx_1);
Result_all.StartBoxON(idx_1) = Result.StartBoxON(idx_1);
Result_all.SentenceON(idx_1) = Result.SentenceON(idx_1);
Result_all.ITION(idx_1) = Result.ITION(idx_1);
Result_all.ProbeON(idx_1) = Result.ProbeON(idx_1);
Result_all.RT(idx_1) = Result.RT(idx_1);
Result_all.CORR(idx_1) = Result.CORR(idx_1);
Result_all.KeyPress(idx_1) = Result.KeyPress(idx_1);
Result_all.WordLocation(idx_1,:,:) = Result.WordLocation(idx_1,:,:);

%%% save out
Para = Para_all;
Result = Result_all;
save([ptbpath subname '.mat'],'cfg','Para','Result')


%% %%%===== connecting the eye data
% just clear the unwanted trial clearly in the end of the first part and begining of the second part
% then copy the second txt into the first one
% because the trigger timepoint counts from the time when the pc is turned
% on and contunuely counting, so it's okay to just put them together


%% %% connecting the MEG data and marker
dataset = [megpath subname filesep subname(3:8) filesep megfile];
cfg         = [];
cfg.dataset = [dataset '.fif'];
hdr         = ft_read_header(cfg.dataset);
event       = ft_read_event(cfg.dataset);
data        = ft_read_data(cfg.dataset);
Trigger_MEG = [[event(strcmp('Trigger',{event.type})).value]' [event(strcmp('Trigger',{event.type})).sample]'];

%%% check how many triggers in the first dataset and remove the extra data and triggers
all_iti_trig = Trigger_MEG(Trigger_MEG(:,1)==16,2);
ntrl = length(all_iti_trig);
end_tp = all_iti_trig(broke_trlid-1)+500; %ITI duration is 500ms
Trigger_MEG(Trigger_MEG(:,2)>end_tp,:) = [];
data(:,end_tp:end) = [];

%%% get the second dataset
cfg         = [];
cfg.dataset = [dataset '-1.fif'];
data2  = ft_read_data(cfg.dataset);
event2 = ft_read_event(cfg.dataset);
Trigger_MEG2 = [[event2(strcmp('Trigger',{event2.type})).value]' [event2(strcmp('Trigger',{event2.type})).sample]'];
Trigger_MEG2(:,2) = Trigger_MEG2(:,2)+end_tp-1;
data = [data data2];
clear data2
Trigger_MEG = [Trigger_MEG; Trigger_MEG2];

% save out
save([savepath 'hdr'], 'hdr','-v7.3')
save([savepath 'data'], 'data','-v7.3')
save([savepath 'Trigger_MEG'], 'Trigger_MEG','-v7.3')



%%% connect dataset that are broken into two parts accidentally
%%% copy from Lexical\ConnectUnexpectedBrokenData

clear

%%% setting parameters
subname = '20211022_b4bc';
megfile = {'b4bc','b4bc-1','b4bc-2'};
ptbfile = {'20211022_b4bc','20211022_b4bc_1'};
savepath = ['Z:\Semantic\Analyse_data\sv_' subname '\'];
ptbpath = 'Z:\Semantic\RawData\PTB_data\';
megpath = 'Z:\Semantic\RawData\MEG_data\';

%% %===== connecting the ptb data
%%% load the second part ptb data
n = length(ptbfile);
trlid_end = [nan(1,n-1) 277];
for ff = n:-1:1
    load([ptbpath ptbfile{ff} '.mat'])
    trlid_start = find(~isnan(Result.FixationON),1);
    disp(['*** broken trl id: ' num2str(trlid_start) '***']);
    if ff == n % use the first-loaded data struct as all data struct
        Para_all = Para;
        Result_all = Result;
    end
    % copy data segment to the all data struct
    idx = trlid_start:trlid_end(ff);
    Para_all.FixDuration(idx) = Para.FixDuration(idx);
    Result_all.FixationON(idx) = Result.FixationON(idx);
    Result_all.StartBoxON(idx) = Result.StartBoxON(idx);
    Result_all.SentenceON(idx) = Result.SentenceON(idx);
    Result_all.ITION(idx) = Result.ITION(idx);
    Result_all.ProbeON(idx) = Result.ProbeON(idx);
    Result_all.RT(idx) = Result.RT(idx);
    Result_all.CORR(idx) = Result.CORR(idx);
    Result_all.KeyPress(idx) = Result.KeyPress(idx);
    Result_all.WordLocation(idx,:,:) = Result.WordLocation(idx,:,:);
    Result_all.EYEdata(idx_1) = Result.EYEdata(idx_1);

    % get the trlid_end
    if ff ~= 1
        trlid_end(ff-1) = trlid_start-1;
    end
end

%%% save out
Para = Para_all;
Result = Result_all;
save([ptbpath subname '.mat'],'cfg','Para','Result')
delete([ptbpath subname '_*.mat'])

%% connecting the MEG data and marker! MEG will break whole session into 2 parts automatically!
dataset = [megpath subname filesep subname(3:8) filesep];
data = [];
Trigger_MEG = [];
trlid_num = [trlid_end(1) diff(trlid_end)];% number of trials in each MEG data segment
if length(trlid_num) ~= length(megfile)
    error('more meg files than expected, need to change trlid_num')
end
add_tp = [0];%time points that needs to be added to trigger_mat,initialize as 0
for ff = 1:length(megfile)
    cfg         = [];
    cfg.dataset = [dataset megfile{ff} '.fif'];
    event       = ft_read_event(cfg.dataset);
    data_tmp    = ft_read_data(cfg.dataset);
    Trig_tmp    = [[event(strcmp('Trigger',{event.type})).value]' [event(strcmp('Trigger',{event.type})).sample]'];
    %%% check how many triggers in this dataset and remove the extra data
    %%% and triggers based on the trlid_end
    all_iti_trig = Trig_tmp(Trig_tmp(:,1)==16,2);
    end_tp = all_iti_trig(trlid_num(ff))+500; %ITI duration is 500ms
    Trig_tmp(Trig_tmp(:,2)>end_tp,:) = [];
    data_tmp(:,end_tp+1:end) = [];
    % put data segment into the all data struct
    data = [data data_tmp];
    Trig_tmp(:,2) = Trig_tmp(:,2)+sum(add_tp);
    Trigger_MEG = [Trigger_MEG; Trig_tmp];
    clear data_tmp Trig_tmp
    % get the to be addded tp from end_tp
    add_tp = [add_tp end_tp];
    % only get hdr for the last segment
    if ff == length(megfile)
        hdr = ft_read_header(cfg.dataset);
    end
end
if length(Trigger_MEG) ~= 1385
    error('wrong number of the triggers')
end
% save out
save([savepath 'hdr'], 'hdr','-v7.3')
save([savepath 'data'], 'data','-v7.3')
save([savepath 'Trigger_MEG'], 'Trigger_MEG','-v7.3')

%% %%%===== connecting the eye data
% after converting the eye-link file from .edf to .asc
% just clear the unwanted trial clearly in the end of the first part and begining of the second part
% then copy the second txt into the first one
% because the trigger timepoint counts from the time when the pc is turned
% on and contunuely counting, so it's okay to just put them together

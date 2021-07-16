% YPan, 20210610
% checking sentences for two projects:
% 1. Semantic violation
% semantic integration between foveal and parafoveal words
% 2. full dynamics of attention during nature reading (from: Degno 2019, JEP)
% using gaze-contingent rapid frequency tagging
% study lexical parafoveal processing + foveal load effect + lingering back
% (spillover effect)

% YPan, 20210621
% making sentences for predictability and plausibility pre-tests


%% semantic violation
fls = 'U:\Projects\Semantic_Parafoveal_Reading\Semantic_violation.txt';
fid = fopen(fls);
SV_Sentmat = {};
Target_ID = []; % word index of targets in each sentence
Words_N = []; % number of words in each sentence
sid = 0; %% index of sentence
SV_PreT_Targ = {}; % get all pre-targets and targets
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        % find all the spaces that seperate words
        space_ids = strfind(strline,' ');
        space_ids = [0 space_ids length(strline)+1];
        targ_id = strfind(strline,'§');
        for wi = 1:length(space_ids)-1
            if space_ids(wi)+1 == targ_id
                istarg = 1;
                Target_ID(sid,1) = wi;
                space_ids(wi) = space_ids(wi)+1;
            else
                istarg = 0;
            end
            SV_Sentmat{sid,wi} = strline(space_ids(wi)+1:space_ids(wi+1)-1);
            %%% get pre-target and target
            if istarg
               SV_PreT_Targ{sid,1} =  SV_Sentmat{sid,wi-1};
               SV_PreT_Targ{sid,2} =  SV_Sentmat{sid,wi};
            end
        end
        Words_N(sid,1) = wi;
    elseif strline == -1
        break
    end
end
fclose(fid);

%%% check whether there're repeated target words
rep_rep = cellfun(@(x) find(strcmp(x,SV_PreT_Targ(:,2))),SV_PreT_Targ(:,2),'Uni',false);
rep_rep = cellfun(@length, rep_rep);
find(rep_rep ~= 1)

%%% 2. get sentence_id with low freq target (CELEX < 10)
% add one emty line then copy into N-watch
SV_PreT_Targ(end+1,:) = {'',''};
% bad target:
!SV_Freq_Leng_Targ = []; %copy the CELEX and Length back here
lowfreq_id = SV_Freq_Leng_Targ(:,1)<10;
badleng_id = SV_Freq_Leng_Targ(:,2)<4 | SV_Freq_Leng_Targ(:,2)>7;
% bad target id in total
bad_targ_id = (lowfreq_id | badleng_id);
length(find(bad_targ_id))

% bad pretarget: copy the CELEX and Length back
!SV_Freq_Leng_PreTarg = [];
lowfreq_id = SV_Freq_Leng_PreTarg(:,1)<10;
badleng_id = SV_Freq_Leng_PreTarg(:,2)<4 | SV_Freq_Leng_PreTarg(:,2)>7;
bad_pretag_id = (lowfreq_id | badleng_id);
length(find(bad_pretag_id))

% total length of pre-target + target
total_leng = SV_Freq_Leng_PreTarg(:,2)+SV_Freq_Leng_Targ(:,2);
figure('name','word length and frequency','color',[1 1 1]);
subplot(2,3,1); histogram(SV_Freq_Leng_PreTarg(:,2),5); title('Pretarget length')
hold on;
subplot(2,3,2); histogram(SV_Freq_Leng_Targ(:,2),4); title('Target length')
subplot(2,3,3); histogram(total_leng,7); title('total length')
subplot(2,1,2); histogram(SV_Freq_Leng_Targ(:,1),50); 
title('InTotal 160 sentences;                    Target frequency')


%%% over all
SV_bad_id = (rep_targ_id | bad_targ_id | bad_pretag_id | tmp);
length(find(SV_bad_id))

%%% remove sentences
find(SV_bad_id)

%% for version_1: first half sentences are incong, the second half are cong
n = size(SV_Sentmat,1);
incong_n = n/2;
Incong_Sent_id = [1:incong_n]';
Cong_Sent_id = [incong_n+1:n]';
% shuffling targets within each pair
tmp_id = [2:2:n/2]';
tmp_id = 2.*tmp_id-1;
tmp_id = kron(tmp_id,[1;1]);
Incong_Sent_id = tmp_id-Incong_Sent_id;
% get new sentences with swapped targets 
new_id = [Incong_Sent_id;Cong_Sent_id];
target_new = SV_PreT_Targ(new_id,2);
SV_Sentmat_V1 = SV_Sentmat;
for i = 1:n
    SV_Sentmat_V1(i,Target_ID(i)) = target_new(i);
end

%% for version_2: first half sentences are cong, the second half are incong
Cong_Sent_id = [1:incong_n]';
Incong_Sent_id = [incong_n+1:n]';
% swapping targets within each pair
tmp_id = [n/2+2:2:n]';
tmp_id = 2.*tmp_id-1;
tmp_id = kron(tmp_id,[1;1]);
Incong_Sent_id = tmp_id-Incong_Sent_id;
% get new sentences with swapped targets 
new_id = [Cong_Sent_id;Incong_Sent_id];
target_new = SV_PreT_Targ(new_id,2);
SV_Sentmat_V2 = SV_Sentmat;
for i = 1:n
    SV_Sentmat_V2(i,Target_ID(i)) = target_new(i);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% make the whole sentence set: SV + WF %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% word frequency sentences from Degno 2019 JEP
fls = 'U:\Projects\Full_attention_dynamics_RIFT\JEP_100sentences.txt';
fid = fopen(fls);
WF_Sentmat = {};
WF_Target_ID = []; % word index of targets in each sentence
WF_Words_N = []; % number of words in each sentence
sid = 0; %% index of sentence
while 1
    strline = fgetl(fid);
    if ischar(strline)
        sid = sid + 1;
        % find all the spaces that seperate words
        space_ids = strfind(strline,' ');
        space_ids = [0 space_ids length(strline)+1];
        targ_id = strfind(strline,'§');
        % only select the starting §
        targ_id = targ_id([1 3]);
        % tmp target_id for each sentence
        tmp_targid = [];
        for wi = 1:length(space_ids)-1
            if ismember(space_ids(wi)+1,targ_id)
                istarg = 1;
                tmp_targid = [tmp_targid wi];
            else
                istarg = 0;
            end
            WF_Sentmat{sid,wi} = strline(space_ids(wi)+1:space_ids(wi+1)-1);
        end
        % get two target_id
        WF_Target_ID(sid,:) = tmp_targid;
        WF_Words_N(sid,1) = wi;
    else
        break
    end
end
fclose(fid);

%%% seperate low and high condition
WF_PreTarg = {};
WF_TargLow = {};
WF_TargHigh = {};
WF_PosTarg = {};
for sid = 1:size(WF_Target_ID,1)
    WF_PreTarg(sid,:) = WF_Sentmat(sid,WF_Target_ID(sid,:)-1);
    WF_PosTarg(sid,:) = WF_Sentmat(sid,WF_Target_ID(sid,:)+1);
    targ = cell2mat(WF_Sentmat(sid,WF_Target_ID(sid,:)));
    targ_id = strfind(targ,'§');
    WF_TargHigh(sid,:) = [{targ(targ_id(1)+1:targ_id(2)-1)},{targ(targ_id(3)+1:targ_id(4)-1)}];
    WF_TargLow(sid,:) = [{targ(targ_id(2)+1:targ_id(3)-1)},{targ(targ_id(4)+1:end)}];
end

%%% V1: the first half low freq, the second half high freq
WF_Sentmat_V1 = WF_Sentmat;
%%% V2: the first half high freq, the second half low freq
WF_Sentmat_V2 = WF_Sentmat;
n = size(WF_Sentmat,1);
for sid = 1:n
    if ~[sid > n/2] % the first half
        WF_Sentmat_V1(sid,WF_Target_ID(sid,:)) = WF_TargLow(sid,:);
        WF_Sentmat_V2(sid,WF_Target_ID(sid,:)) = WF_TargHigh(sid,:);
    else
        WF_Sentmat_V2(sid,WF_Target_ID(sid,:)) = WF_TargLow(sid,:);
        WF_Sentmat_V1(sid,WF_Target_ID(sid,:)) = WF_TargHigh(sid,:);
    end
end

%% SV + WF
% make 2 sentence-sets the same width
SV_Sentmat_V1(:,[16 17]) = {[]};
SV_Sentmat_V2(:,[16 17]) = {[]};
Target_ID(:,2) = nan;
% combine
All_Sentmat_V1 = [SV_Sentmat_V1;WF_Sentmat_V1]; %[incog;cong][low;high]
All_Sentmat_V2 = [SV_Sentmat_V2;WF_Sentmat_V2]; %[cog;incong][high;low]
targid_all = [Target_ID; WF_Target_ID];
% condition id: SV--1(incong--1); WF--2(low---1)
cond_id_V1 = [11.*ones(length(SV_Sentmat_V1)/2,1);12.*ones(length(SV_Sentmat_V1)/2,1);...
              21.*ones(length(WF_Sentmat_V1)/2,1);22.*ones(length(WF_Sentmat_V1)/2,1)];
cond_id_V2 = [12.*ones(length(SV_Sentmat_V1)/2,1);11.*ones(length(SV_Sentmat_V1)/2,1);...
              22.*ones(length(WF_Sentmat_V1)/2,1);21.*ones(length(WF_Sentmat_V1)/2,1)];
% randeromazation
randid = randperm(length(All_Sentmat_V1))';
SentMat_1 = All_Sentmat_V1(randid,:);
SentMat_2 = All_Sentmat_V2(randid,:);
Cond_V1 = cond_id_V1(randid);
Cond_V2 = cond_id_V2(randid);
Rand_Target_ID = targid_all(randid,:);

% get questions
filnam = 'U:\Projects\Semantic_Parafoveal_Reading\questions.xlsx';
[~,txt] = xlsread(filnam,'S1:T260'); % random order sentmat
% remove the start and end cell from txt (rows of txt might be smaller than 
% sentences number, we need to pre-occupy start and end cell)
txt(1,1) = txt(1,2);
txt(end,2) = txt(end,1);
Question = txt;

% check questions
id = ~cellfun(@isempty,Question(:,1));
s1 = SentMat(id,:);
q = Question(id,:);

% saveout
SentMat_all.FixOrder_V1 = All_Sentmat_V1;
SentMat_all.FixOrder_V2 = All_Sentmat_V2;
SentMat_all.Words_N = [Words_N;WF_Words_N];
SentMat_all.RandId = randid;
SentMat_all.Rand_CondID_hdr = {'11-SV_Incong','12-SV_Cong','21-WF_Low','22-WF_High'};
SentMat_all.Rand_CondID_V1 = Cond_V1;
SentMat_all.Rand_CondID_V2 = Cond_V2;
SentMat_all.Rand_Target_ID = Rand_Target_ID;
SentMat_all.Rand_Question = Question;
save('U:\Projects\Semantic_Parafoveal_Reading\SentMat_all.mat','SentMat_all');


% for PTB, needs: Sentmat_V1/2, Cond_V1/2, Rand_Target_ID
cd Z:\Semantic\PTB_codes\
SentMat = SentMat_1;
save SentMat_1 SentMat
SentMat = SentMat_2;
save SentMat_2 SentMat
save Question Question
TargLoc_Cond = [Rand_Target_ID Cond_V1];
save TargLoc_Cond_V1 TargLoc_Cond
TargLoc_Cond = [Rand_Target_ID Cond_V2];
save TargLoc_Cond_V2 TargLoc_Cond



%% 20210716: select word_freq sentences out of the semantiv task, put orth_freq sentences in as fillers
load('SentMat_2.mat')
load('TargLoc_Cond_2.mat')
sv_2_id = TargLoc_Cond(:,3)<20;
sentmat_sv_2 = SentMat(sv_2_id,:);
wf_2_id = TargLoc_Cond(:,3)>20;
sentmat_wf_2 = SentMat(wf_2_id,:);
cond_sv_2 = TargLoc_Cond(sv_2_id,:);
cond_wf_2 = TargLoc_Cond(wf_2_id,:);
Question_sv = Question(sv_2_id,:);
Question_wf = Question(wf_2_id,:);

%%
SentMat = sentmat_wf_1;
TargLoc_Cond = cond_wf_1;
save FA_SentMat_1 SentMat
save FA_TargLoc_Cond_1 TargLoc_Cond
%%
SentMat = sentmat_wf_2;
TargLoc_Cond = cond_wf_2;
save FA_SentMat_2 SentMat
save FA_TargLoc_Cond_2 TargLoc_Cond
%% 
Question = Question_wf;
save FA_Question Question

%
allid = allid == 1;
cond_sv_1 = cond_sv_1(:,[1 3 2]);
cond_sv_2 = cond_sv_2(:,[1 3 2]);

%%
allsent_1 = cell(277,17);
cond_1 = zeros(277,3);
question_1 = cell(277,2);

allsent_1(allid,:) = sentmat_sv_1;
cond_1(allid,:) = cond_sv_1;
question_1(allid,:) = Question_sv;

allsent_1(~allid,1:16) = sentmat_of_1;
cond_1(~allid,:) = TargCondPair_of_1;
question_1(~allid,:) = Question_of_1;

%
allsent_2 = cell(277,17);
cond_2 = zeros(277,3);
question_2 = cell(277,2);

allsent_2(allid,:) = sentmat_sv_1;
cond_2(allid,:) = cond_sv_1;
question_2(allid,:) = Question_sv;

allsent_2(~allid,1:16) = sentmat_of_2;
cond_2(~allid,:) = TargCondPair_of_2;
question_2(~allid,:) = Question_of_2;

%
allsent_3 = cell(277,17);
cond_3 = zeros(277,3);
question_3 = cell(277,2);

allsent_3(allid,:) = sentmat_sv_1;
cond_3(allid,:) = cond_sv_1;
question_3(allid,:) = Question_sv;

allsent_3(~allid,1:16) = sentmat_of_3;
cond_3(~allid,:) = TargCondPair_of_3;
question_3(~allid,:) = Question_of_3;

%
allsent_4 = cell(277,17);
cond_4 = zeros(277,3);
question_4 = cell(277,2);

allsent_4(allid,:) = sentmat_sv_2;
cond_4(allid,:) = cond_sv_2;
question_4(allid,:) = Question_sv;

allsent_4(~allid,1:16) = sentmat_of_1;
cond_4(~allid,:) = TargCondPair_of_1;
question_4(~allid,:) = Question_of_1;

%
allsent_5 = cell(277,17);
cond_5 = zeros(277,3);
question_5 = cell(277,2);

allsent_5(allid,:) = sentmat_sv_2;
cond_5(allid,:) = cond_sv_2;
question_5(allid,:) = Question_sv;

allsent_5(~allid,1:16) = sentmat_of_2;
cond_5(~allid,:) = TargCondPair_of_2;
question_5(~allid,:) = Question_of_2;

%
allsent_6 = cell(277,17);
cond_6 = zeros(277,3);
question_6 = cell(277,2);

allsent_6(allid,:) = sentmat_sv_2;
cond_6(allid,:) = cond_sv_2;
question_6(allid,:) = Question_sv;

allsent_6(~allid,1:16) = sentmat_of_3;
cond_6(~allid,:) = TargCondPair_of_3;
question_6(~allid,:) = Question_of_3;

%%
for i = 1:6
    eval(['SentMat = allsent_' num2str(i) ';']);
    eval(['TargCondPair = cond_' num2str(i) ';']);
    eval(['Question = question_' num2str(i) ';']);
    save(['SentMat_' num2str(i)], 'SentMat')
    save(['TargCondPair_' num2str(i)], 'TargCondPair')
    save(['Question_' num2str(i)], 'Question')
end

































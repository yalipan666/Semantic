% YPan. 20210703
% Making the sentences for pretest: predictability and plausibility
% Analyzing the pretests
% (first do predictability test to make sure no highly predicted targets, 
% then do the plausibility test twice for two versions)

%%% path
ppath = 'U:\Projects\Semantic_Parafoveal_Reading\Pre-test\';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% making pre-tests for semantic violation %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prediction test
fls = 'U:\Projects\Semantic_Parafoveal_Reading\Semantic_violation.txt';
fid = fopen(fls);
sid = 0; %% index of sentence
Sent_before_targ = {};  
Target = {}; 
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        % find target words
        targ_id = strfind(strline,'§');
        space_ids = strfind(strline,' ');
        space_ids = space_ids(space_ids > targ_id);
        Target{sid,1} = strline(targ_id+1:space_ids(1)-1);
        Sent_before_targ{sid,1} = [strline(1:targ_id-1) '..................'];
    elseif strline == -1
        break
    end
end
fclose(fid);

% randermization
randid = randperm(size(Sent_before_targ,1));
pred_sents = Sent_before_targ(randid);
pred_targers = Target(randid);

save randid randid

%% plausibility test
% get all the congruent sentences
fls = 'U:\Projects\Semantic_Parafoveal_Reading\Semantic_violation.txt';
fid = fopen(fls);
SV_Sentmat = {};
SV_PreTarg_Targ_PosTarg = {};
Target_ID = []; % word index of targets in each sentence
Words_N = []; % number of words in each sentence
sid = 0; %% index of sentence
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        %%% find all the spaces that seperate words
        % remove the ocassional empty space at the end of sentences
        len_strline = strfind(strline,'.');
        strline = strline(1:len_strline);
        space_ids = strfind(strline,' ');
        space_ids = [0 space_ids len_strline+1];
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
        end
        Words_N(sid,1) = wi;
        %%% get pre-target, target, and post-target
        SV_PreTarg_Targ_PosTarg{sid,1} =  SV_Sentmat{sid,Target_ID(sid,1)-1};
        SV_PreTarg_Targ_PosTarg{sid,2} =  SV_Sentmat{sid,Target_ID(sid,1)};
        SV_PreTarg_Targ_PosTarg{sid,3} =  SV_Sentmat{sid,Target_ID(sid,1)+1};
        
    elseif strline == -1
        break
    end
end
fclose(fid);

%%% shuffle targets to create incongruent sentences for the first half of
%%% sentences
n = size(SV_Sentmat,1);
incong_n = n/2;

%% for version_1: first half sentences are incong, the second half are cong
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
SV_Sentmat_new = SV_Sentmat;
for i = 1:n
    SV_Sentmat_new(i,Target_ID(i)) = target_new(i);
end

%%% load in sentences with middle plausibility
fls = [ppath 'plausibility\middle_plausible.txt'];
fid = fopen(fls);
sid = 0; %% index of sentence
Sent_MidPlaus = {};  
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        Sent_MidPlaus{sid,1} = strline;
    elseif strline == -1
        break
    end
end
fclose(fid);
MidPlaus_Sent_id = n+[1:length(Sent_MidPlaus)]';
% make setences in SV_Sentmat_new into oneline string
SV_Sentmat_string = {};
for i = 1:n
    tmp_string = [];
    for w = 1:length(SV_Sentmat_new(i,:))
        tmp_string = [tmp_string SV_Sentmat_new{i,w} ' '];
    end
    SV_Sentmat_string{i,1} = tmp_string;
end 
All_strings = [SV_Sentmat_string; Sent_MidPlaus];
% randermization
all_sent_id = [Incong_Sent_id; Cong_Sent_id; MidPlaus_Sent_id];
randid = randperm(length(all_sent_id))';
% randid = zeros(230,1);
% randid(PreTest.Plausibility.incongID_v2) = Incong_Sent_id(randperm(80));
% randid(PreTest.Plausibility.congID_v2) = Cong_Sent_id(randperm(80));
% mid_loc = setdiff(1:230,[PreTest.Plausibility.incongID_v2;PreTest.Plausibility.congID_v2]);
% mid_id = 161:230;
% randid(mid_loc) = mid_id(randperm(70));

% check whether more than 3 sentences from the condition are in a row
sv_id = randid<n+1;
incong_id = randid<incong_n+1;
cong_id = randid>incong_n & randid<n+1;
Rand_strings = All_strings(randid);
%%% then copy Rand_strings to a word document, remove all the '' at the
%%% begining and end of sentences (please not cases like ['s,don't,'d] where ' is
%%% needed). Finally, copy sentences to the excel.
PreTest.Plausibility.incongID_v1 = find(incong_id);
PreTest.Plausibility.congID_v1 = find(cong_id);

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
SV_Sentmat_new = SV_Sentmat;
for i = 1:n
    SV_Sentmat_new(i,Target_ID(i)) = target_new(i);
end

%%% load in sentences with middle plausibility
fls = [ppath 'plausibility\middle_plausible.txt'];
fid = fopen(fls);
sid = 0; %% index of sentence
Sent_MidPlaus = {};  
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        Sent_MidPlaus{sid,1} = strline;
    elseif strline == -1
        break
    end
end
fclose(fid);
MidPlaus_Sent_id = n+[1:length(Sent_MidPlaus)]';
% make setences in SV_Sentmat_new into oneline string
SV_Sentmat_string = {};
for i = 1:n
    tmp_string = [];
    for w = 1:length(SV_Sentmat_new(i,:))
        tmp_string = [tmp_string SV_Sentmat_new{i,w} ' '];
    end
    SV_Sentmat_string{i,1} = tmp_string;
end 
All_strings = [SV_Sentmat_string; Sent_MidPlaus];
% randermization
all_sent_id = [Incong_Sent_id; Cong_Sent_id; MidPlaus_Sent_id];
% randid = randperm(length(all_sent_id))';
randid = zeros(230,1);
randid(PreTest.Plausibility.incongID_v2) = Incong_Sent_id(randperm(80));
randid(PreTest.Plausibility.congID_v2) = Cong_Sent_id(randperm(80));
mid_loc = setdiff(1:230,[PreTest.Plausibility.incongID_v2;PreTest.Plausibility.congID_v2]);
mid_id = 161:230;
randid(mid_loc) = mid_id(randperm(70));

% check whether more than 3 sentences from the condition are in a row
sv_id = randid<n+1;
incong_id = randid>incong_n & randid<n+1;
cong_id = randid<incong_n+1;
Rand_strings = All_strings(randid);
%%% then copy Rand_strings to a word document, remove all the '' at the
%%% begining and end of sentences (please not cases like ['s,don't,'d], where ' is
%%% needed). Finally, copy sentences to the excel.
PreTest.Plausibility.incongID_v2 = find(incong_id);
PreTest.Plausibility.congID_v2 = find(cong_id);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% analyzing pre-tests for semantic violation %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1st prediction test for the origianl 158 sentences
path_pred = [ppath 'predictability\subject_files\first_158sents\'];
%%% load in sentences of predictability test
fls = [path_pred 'Predictability_test_NoQuoteMark.txt'];
fid = fopen(fls);
sent_nodot = {};  
sid = 0;
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        % delete ..... from sentences
        strline = erase(strline,'..................');
        sent_nodot{sid,1} = strline;
    elseif strline == -1
        break
    end
end
fclose(fid);

%!!!!manually copy sentences into one big .txt 
%%% load in sentences Predictability_all_in_one
fls = [path_pred 'Predictability_allsubs_NoQuoteMark.txt'];
fid = fopen(fls);
n = length(sent_nodot);
predict_targ = cell(n,1);  
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        % find sentence_id 
        dot_ids = strfind(strline,'.');
        % since the id contiues across subs, we need get mod as the real id
        sid = str2double(strline(1:dot_ids(1)-1));
        % get the first word from answers
        strline = erase(strline(dot_ids(1)+2:end),sent_nodot{sid,1});
        strline = erase(strline,'.');
        % remove the space between answer by iteration
        while length(strline)>1
            if strcmp(strline(1),' ')
                strline = strline(2:end);
            else
                break
            end
        end
        % get the first word from answer
        space_ids = strfind(strline,' ');
        if ~isempty(space_ids) % more than one word in the answer
            strline = strline(1:space_ids-1);
        end
        predict_targ{sid,1} = [predict_targ{sid,1};{strline}];
    elseif strline == -1
        break
    end
end
fclose(fid);
% check if there's empty cell
% get the prediction rate of the real target word
PreTest.Predictability = cell(n+1,6);
PreTest.Predictability(1,:) = {'Sent4test','target','percent4target',...
    'most_predict_word','percent4mosprewrd','Nsub4thisItem'};
PreTest.Predictability(2:end,[1 2]) = [sent_nodot pred_targers];
for i = 1:n
    ans_wrd = tabulate(predict_targ{i, 1});  
    targid = find(strcmp(pred_targers{i,1},ans_wrd(:,1)));
    if isempty(targid)
        PreTest.Predictability{i+1,3} = 0;
    else
        PreTest.Predictability{i+1,3} = ans_wrd{targid,3};
    end
    [~,most_id] = max(cell2mat(ans_wrd(:,2)));
    PreTest.Predictability(i+1,[4 5]) = ans_wrd(most_id,[1 3]);
    PreTest.Predictability{i+1,6} = sum(cell2mat(ans_wrd(:,2)));
end
% copy out as csv --Predictability_result(sheet:1st) (csvwrite doesn't accept cell array)
 
% (details of the revision of the original sentences, please see changes from 1st pretest.doc )

%% 2nd prediction test for the revised 8 sentences
path_pred = [ppath 'predictability\subject_files\second_8sents\'];
%%% load in sentences of predictability test
fls = [path_pred 'Predictability_test_NoQuoteMark_2.txt'];
fid = fopen(fls);
sent_nodot = {};  
sid = 0;
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        % delete ..... from sentences
        strline = erase(strline,'..................');
        sent_nodot{sid,1} = strline;
    elseif strline == -1
        break
    end
end
fclose(fid);

%!!!!manually copy sentences into one big .txt 
%%% load in sentences Predictability_all_in_one
fls = [path_pred 'Predictability_allsubs_NoQuoteMark_2.txt'];
fid = fopen(fls);
n = length(sent_nodot);
predict_targ = cell(n,1);  
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        % find sentence_id 
        dot_ids = strfind(strline,'.');
        % since the id contiues across subs, we need get mod as the real id
        sid = str2double(strline(1:dot_ids(1)-1));
        % get the first word from answers
        strline = erase(strline(dot_ids(1)+2:end),sent_nodot{sid,1});
        strline = erase(strline,'.');
        % remove the space between answer by iteration
        while length(strline)>1
            if strcmp(strline(1),' ')
                strline = strline(2:end);
            else
                break
            end
        end
        % get the first word from answer
        space_ids = strfind(strline,' ');
        if ~isempty(space_ids) % more than one word in the answer
            strline = strline(1:space_ids-1);
        end
        predict_targ{sid,1} = [predict_targ{sid,1};{strline}];
    elseif strline == -1
        break
    end
end
fclose(fid);
% check if there's empty cell
% get the prediction rate of the real target word
pred_targers = {'nurse';'report';'drone';'coffee';'crisis';'cups';'paper';'baker'};
PreTest.Predictability_2 = cell(n+1,6);
PreTest.Predictability_2(1,:) = {'Sent4test','target','percent4target',...
    'most_predict_word','percent4mosprewrd','Nsub4thisItem'};
PreTest.Predictability_2(2:end,[1 2]) = [sent_nodot pred_targers];
for i = 1:n
    ans_wrd = tabulate(predict_targ{i, 1});  
    targid = find(strcmp(pred_targers{i,1},ans_wrd(:,1)));
    if isempty(targid)
        PreTest.Predictability_2{i+1,3} = 0;
    else
        PreTest.Predictability_2{i+1,3} = ans_wrd{targid,3};
    end
    [~,most_id] = max(cell2mat(ans_wrd(:,2)));
    PreTest.Predictability_2(i+1,[4 5]) = ans_wrd(most_id,[1 3]);
    PreTest.Predictability_2{i+1,6} = sum(cell2mat(ans_wrd(:,2)));
end
% copy out as csv --Predictability_result(sheet:2nd) (csvwrite doesn't accept cell array)


%% ============ plausibility test ============= %%
%%%% V1--first_old158sents
% get all files
plau_path = [ppath 'plausibility\subject_files\first_old158sents\'];
MyFolderInfo = dir(plau_path);
n_sub = length(MyFolderInfo)-2;
plaus_rate = zeros(228,n_sub);
for i = 1:n_sub
    filnam = MyFolderInfo(i+2).name;
    ratings = xlsread([plau_path filnam],1,'G16:G255'); 
    plaus_rate(:,i) = ratings(~isnan(ratings));
end
% seperate conditions
incong_id = PreTest.Plausibility.incongID_v1_old158;
cong_id = PreTest.Plausibility.congID_v1_old158;
% get mean and std over subs
all_cong = plaus_rate(cong_id,:);
sub_cong = mean(all_cong,1);
all_incong = plaus_rate(incong_id,:);
sub_incong = mean(all_incong,1);
% check plausibility item-wise
% get whole sentence mat
[~,txt] = xlsread([plau_path filnam],1,'A16:A255'); 
txt = txt(~cellfun(@isempty,txt));
item_cong = mean(all_cong,2);
low_cong_id = cong_id(item_cong < 4)
item_cong(item_cong < 4)
txt(low_cong_id)
item_incong = mean(all_incong,2);
low_incong_id = incong_id(item_incong > 3)
item_cong(item_incong > 3)
txt(low_incong_id)

% initializing for the combined data for version1
sub_cong_all = sub_cong; % 1*n 
sub_incong_all = sub_incong;

%%%% V1--first_old160sents
% get all files
plau_path = [ppath 'plausibility\subject_files\first_old160sents\'];
MyFolderInfo = dir(plau_path);
n_sub = length(MyFolderInfo)-2;
plaus_rate = zeros(230,n_sub);
for i = 1:n_sub
    filnam = MyFolderInfo(i+2).name;
    ratings = xlsread([plau_path filnam],1,'G16:G257'); 
    plaus_rate(:,i) = ratings(~isnan(ratings));
end
% seperate conditions
incong_id = PreTest.Plausibility.incongID_v1_old160;
cong_id = PreTest.Plausibility.congID_v1_old160;
% get mean and std over subs
all_cong = plaus_rate(cong_id,:);
sub_cong = mean(all_cong,1);
all_incong = plaus_rate(incong_id,:);
sub_incong = mean(all_incong,1);
% check plausibility item-wise
% get whole sentence mat
[~,txt] = xlsread([plau_path filnam],1,'A16:A257'); 
txt = txt(~cellfun(@isempty,txt));
item_cong = mean(all_cong,2);
low_cong_id = cong_id(item_cong < 4)
item_cong(item_cong < 4)
txt(low_cong_id)
item_incong = mean(all_incong,2);
low_incong_id = incong_id(item_incong > 3)
item_cong(item_incong > 3)
txt(incong_id(low_incong_id))

%%% combine together --- first_old158sents & first_old160sents
sub_cong_all = [sub_cong_all sub_cong]; % 1*n 
sub_incong_all = [sub_incong_all sub_incong];


%% === new-version-1
plau_path = [ppath 'plausibility\subject_files\version_1\'];
MyFolderInfo = dir(plau_path);
n_sub = length(MyFolderInfo)-2;
plaus_rate = zeros(230,n_sub);
for i = 1:n_sub
    filnam = MyFolderInfo(i+2).name;
    ratings = xlsread([plau_path filnam],1,'G16:G255'); 
    plaus_rate(:,i) = ratings(~isnan(ratings));
end
% seperate conditions
incong_id = PreTest.Plausibility.incongID_v1;
cong_id = PreTest.Plausibility.congID_v1;
% get mean and std over subs
all_cong = plaus_rate(cong_id,:);
sub_cong = mean(all_cong,1);
all_incong = plaus_rate(incong_id,:);
sub_incong = mean(all_incong,1);
% check plausibility item-wise
% get whole sentence mat
[~,txt] = xlsread([plau_path filnam],1,'A16:A255'); 
txt = txt(~cellfun(@isempty,txt));
item_cong = mean(all_cong,2);
low_cong_id = cong_id(item_cong < 4)
item_cong(item_cong < 4)
txt(low_cong_id) 
item_incong = mean(all_incong,2);
low_incong_id = incong_id(item_incong > 3)
item_incong(item_incong > 3)
txt(low_incong_id)

%%% combine together --- V1
sub_cong_v1 = [sub_cong_all sub_cong]; % 1*n 
sub_incong_v1 = [sub_incong_all sub_incong];
sub_cong = sub_cong_v1;
sub_incong = sub_incong_v1;

%%%% save out to pre-test struct
PreTest.PlausibilityResult.V1_age_female = ...
    [22 26 23 23 22 22 23 23 19 31 23 25 50 22 23 23 23 19 19 20 18 19 19 19 19 20 20;
     0 0 nan 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1];
PreTest.PlausibilityResult.V1_cong_subs = sub_cong;
PreTest.PlausibilityResult.V1_cong_avg_std = [mean(sub_cong) std(sub_cong)]; 
PreTest.PlausibilityResult.V1_incong_subs = sub_incong;
PreTest.PlausibilityResult.V1_incong_avg_std = [mean(sub_incong) std(sub_incong)];  
save([ppath 'PreTest.mat'],'PreTest');



%% === new-version-2
ppath = 'U:\Projects\Semantic_Parafoveal_Reading\Pre-test\';
load('U:\Projects\Semantic_Parafoveal_Reading\Pre-test\PreTest');
plau_path = [ppath 'plausibility\subject_files\version_2\'];
MyFolderInfo = dir(plau_path);
n_sub = length(MyFolderInfo)-2;
plaus_rate = zeros(230,n_sub);
for i = 1:n_sub
    filnam = MyFolderInfo(i+2).name;
    ratings = xlsread([plau_path filnam],1,'G16:G255'); 
    plaus_rate(:,i) = ratings(~isnan(ratings));
end
% seperate conditions
incong_id = PreTest.Plausibility.incongID_v2;
cong_id = PreTest.Plausibility.congID_v2;
% get mean and std over subs
all_cong = plaus_rate(cong_id,:);
sub_cong = mean(all_cong,1);
all_incong = plaus_rate(incong_id,:);
sub_incong = mean(all_incong,1);
% save out
PreTest.PlausibilityResult.V2_age_female = ...
    [26 25 23 23 24 24 21 23 19 20 20 20 21 20 20 19 19 19 18 19 20;
     0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
PreTest.PlausibilityResult.V2_cong_subs = sub_cong;
PreTest.PlausibilityResult.V2_cong_avg_std = [mean(sub_cong) std(sub_cong)]; 
PreTest.PlausibilityResult.V2_incong_subs = sub_incong;
PreTest.PlausibilityResult.V2_incong_avg_std = [mean(sub_incong) std(sub_incong)];  
save([ppath 'PreTest.mat'],'PreTest');


% check plausibility item-wise
% get whole sentence mat
[~,txt] = xlsread([plau_path filnam],1,'A16:A255'); 
txt = txt(~cellfun(@isempty,txt));
item_cong = mean(all_cong,2);
low_cong_id = cong_id(item_cong < 4)
item_cong(item_cong < 4)
txt(low_cong_id) % bad sentences
item_incong = mean(all_incong,2);
low_incong_id = incong_id(item_incong > 3)
item_incong(item_incong > 3)
txt(low_incong_id)

sub_cong_v2 = sub_cong; % 1*n 
sub_incong_v2 = sub_incong;


%%
[mean(sub_cong_v1) std(sub_cong_v1)]
[mean(sub_incong_v1) std(sub_incong_v1)]

[mean(sub_cong_v2) std(sub_cong_v2)]
[mean(sub_incong_v2) std(sub_incong_v2)]


[mean(sub_cong_v1([1:10 18:24])) std(sub_cong_v1([1:10 18:24]))]
[mean(sub_incong_v1([1:10 18:24])) std(sub_incong_v1([1:10 18:24]))]


PreTest.Plausibility.incong_avg_std = [mean(all_incong_rate) std(all_incong_rate)];  
PreTest.Plausibility.cong_avg_std = [mean(all_cong_rate) std(all_cong_rate)];  
save([ppath 'PreTest.mat'],'PreTest');


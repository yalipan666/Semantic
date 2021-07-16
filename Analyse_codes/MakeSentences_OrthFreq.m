%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% make sentences for OrthWrdFrq project %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sentences are from 
% White, S. J. (2008). Eye movement control during reading: Effects of word
% frequency and orthographic familiarity. Journal of experimental psychology: 
% Human perception and performance, 34(1), 205.
% Sentences a, b, and c refer to the frequent and orthographically familiar,
% infrequent and orthographically familiar, and infrequent and orthographically
% unfamiliar conditions, respectively.
ppath = 'U:\Projects\Orth_Freq\';
fls = [ppath 'orth_freq.txt'];
fid = fopen(fls);
OF_SentMat = {};
Target_ID = []; % word index of targets in each sentence
Words_N = []; % number of words in each sentence
sid = 0; %% index of sentence
pairid_condid = [];
while 1
    strline = fgetl(fid);
    if ischar(strline) && ~isempty(strline)
        sid = sid + 1;
        % get the condition string
        dot_id = strfind(strline,'.');
        dot_id = dot_id(1);
        pairid_condid(sid,1) = str2double(strline(1:dot_id-2));
        pairid_condid(sid,2) = strfind('abc',strline(dot_id-1));
        % find all the spaces that seperate words
        strline = strline(dot_id+2:end);
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
            OF_SentMat{sid,wi} = strline(space_ids(wi)+1:space_ids(wi+1)-1);
        end
        Words_N(sid,1) = wi;
    elseif strline == -1
        break
    end
end
fclose(fid);

%% version-1
% pseudo-randomize sentences
n_pair = length(unique(pairid_condid(:,1)));
% randomize sentence_pair in every third of the sentences
pair_id = [randperm(n_pair)' ;randperm(n_pair)'; randperm(n_pair)' ];
% randomize conditions
condid_3rd_1 = pairid_condid(1:n_pair,2);
condid_3rd_1 = condid_3rd_1(randperm(n_pair));
% check it to make sure no more than 3 sentences are of the same condition in a row
cond_id = zeros(size(pair_id));
randid = zeros(size(pair_id));
pid = pairid_condid(:,1);
cid = pairid_condid(:,2);
for i = 1:n_pair
    tmp_p = find(pair_id == i);
    tmp_c = randperm(3)';
    c1 = condid_3rd_1(tmp_p(1));
    tmp_c = setdiff(tmp_c,c1,'stable');
    c2 = tmp_c(1);
    c3 = tmp_c(2);
    cond_id(tmp_p) = [c1;c2;c3];
    % get the pseudo-randomize index
    p = pid==i;
    randid(tmp_p) = [find(p&cid==c1);find(p&cid==c2);find(p&cid==c3)];
end
SentMat_rand = OF_SentMat(randid,:);
Target_ID_rand = Target_ID(randid);
condid_pairid_rand = [cond_id pair_id];
Orth_Freq.SentMat_1 = SentMat_rand;
Orth_Freq.Target_ID = Target_ID_rand;
Orth_Freq.condid_hdr = {'1:freq_fami','infreq_fami','infreq_unfami'};
Orth_Freq.condid_pairid = condid_pairid_rand;
Orth_Freq.Words_N = Words_N;

% saveout
SentMat = SentMat_rand;
TargCondPair = [Target_ID_rand condid_pairid_rand];
save('Z:\Orth_Freq\PTB_codes\SentMat_1','SentMat');
save('Z:\Orth_Freq\PTB_codes\TargCondPair_1','TargCondPair');


%% version-2 & 3
randid_2 = zeros(size(randid));
randid_3 = zeros(size(randid));
cond_id_2 = zeros(size(randid));
cond_id_3 = zeros(size(randid));
pid = pairid_condid(:,1); % paid id for raw SentMat
cid = pairid_condid(:,2); % cond id for raw SentMat
allcond = [2 3;3 1;1 2];
for s = 1:size(SentMat_rand,1)
    c1 = condid_pairid_rand(s,1);
    c2 = allcond(c1,1);
    c3 = allcond(c1,2);
    cond_id_2(s) = c2;
    cond_id_3(s) = c3;
    % get the pseudo-randomize index
    p = pid==condid_pairid_rand(s,2);
    randid_2(s) = find(p&cid==c2);
    randid_3(s) = find(p&cid==c3);
end
SentMat_2 = OF_SentMat(randid_2,:);
SentMat_3 = OF_SentMat(randid_3,:);
condid_pairid_2 = [cond_id_2 pair_id];
condid_pairid_3 = [cond_id_3 pair_id];
% % % check whether some cond_pair are repeated
% % tmp = condid_pairid_2(:,1).*100 + condid_pairid_2(:,2);
% % length(unique(tmp))
% % tmp = condid_pairid_3(:,1).*100 + condid_pairid_3(:,2);
% % length(unique(tmp))
Orth_Freq.SentMat_2 = SentMat_2;
Orth_Freq.condid_pairid_2 = condid_pairid_2;
Orth_Freq.SentMat_3 = SentMat_3;
Orth_Freq.condid_pairid_3 = condid_pairid_3;
save([ppath 'Orth_Freq'],'Orth_Freq')

% saveout
SentMat = SentMat_2;
TargCondPair = [Target_ID_rand condid_pairid_2];
save('Z:\Orth_Freq\PTB_codes\SentMat_2','SentMat');
save('Z:\Orth_Freq\PTB_codes\TargCondPair_2','TargCondPair');
SentMat = SentMat_3;
TargCondPair = [Target_ID_rand condid_pairid_3];
save('Z:\Orth_Freq\PTB_codes\SentMat_3','SentMat');
save('Z:\Orth_Freq\PTB_codes\TargCondPair_3','TargCondPair');


%%%%  make questions for each versions
%% version 1: make setences into oneline string
SentMat_string = {};
for i = 1:size(SentMat_rand,1)
    tmp_string = [];
    for w = 1:length(SentMat_rand(i,:))
        tmp_string = [tmp_string SentMat_rand{i,w} ' '];
    end
    SentMat_string{i,1} = tmp_string;
end 
% use these sentences as a reference 
filnam = [ppath 'Questions_1.xlsx'];
[~,txt] = xlsread(filnam,'A1:B117'); % random order SentMat
% remove the start and end cell from txt (rows of txt might be smaller than 
% sentences number, we need to pre-occupy start and end cell)
txt(1,1) = txt(1,2);
txt(end,2) = txt(end,1);
Question = txt;
% check questions
id = ~cellfun(@isempty,Question(:,1));
s1 = SentMat_string(id,:);
q = Question(id,:);
[s1 q];
save('Z:\Orth_Freq\PTB_codes\Question_1','Question');

%% version 2: make setences into oneline string
SentMat_string = {};
for i = 1:size(SentMat_2,1)
    tmp_string = [];
    for w = 1:length(SentMat_2(i,:))
        tmp_string = [tmp_string SentMat_2{i,w} ' '];
    end
    SentMat_string{i,1} = tmp_string;
end 
% use these sentences as a reference 
filnam = [ppath 'Questions_2.xlsx'];
[~,txt] = xlsread(filnam,'A1:B117'); % random order SentMat
% remove the start and end cell from txt (rows of txt might be smaller than 
% sentences number, we need to pre-occupy start and end cell)
txt(1,1) = txt(1,2);
Question = txt;
% check questions
id = ~cellfun(@isempty,Question(:,1));
s1 = SentMat_string(id,:);
q = Question(id,:);
[s1 q];
save('Z:\Orth_Freq\PTB_codes\Question_2','Question');

%% version 3: make setences into oneline string
SentMat_string = {};
for i = 1:size(SentMat_3,1)
    tmp_string = [];
    for w = 1:length(SentMat_3(i,:))
        tmp_string = [tmp_string SentMat_3{i,w} ' '];
    end
    SentMat_string{i,1} = tmp_string;
end 
% use these sentences as a reference 
filnam = [ppath 'Questions_3.xlsx'];
[~,txt] = xlsread(filnam,'A1:B117'); % random order SentMat
% remove the start and end cell from txt (rows of txt might be smaller than 
% sentences number, we need to pre-occupy start and end cell)
txt(1,1) = txt(1,2);
txt(end,2) = txt(end,1);
Question = txt;
% check questions
id = ~cellfun(@isempty,Question(:,1));
s1 = SentMat_string(id,:);
q = Question(id,:);
[s1 q];
save('Z:\Orth_Freq\PTB_codes\Question_3','Question');



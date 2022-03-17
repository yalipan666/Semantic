%% get the statistic information about the word properties of words in sv task

%% get the words
cd('Z:\Semantic\PTB_codes\');
load('SentMat_1.mat')
load('TargCondPair_1.mat')
sv = TargCondPair(:,2)>10; %index for sv
target_loc = TargCondPair(sv,1);
sv_mat = SentMat(sv,:);
n = size(sv_mat,1);
pretarget = cell(n,1);
target = cell(n,1);
postarget = cell(n,1);
for i = 1:n
    pretarget{i,1} = sv_mat{i,target_loc(i)-1};
    target{i,1} = sv_mat{i,target_loc(i)};
    postarget{i,1} = sv_mat{i,target_loc(i)+1};
end
words = [pretarget;target;postarget];
words{n+1} = 'nan'; %add one extra line to the words because the N-watch always ignored the last word when paste the word list


%% copy the list of words into N-watch, select 'CELEXT total' and 'Len_L'

%% copy back the freq and length
freq = [];
len = [];

freq_pre_tar_pos = [freq(1:n)' freq(n+1:2*n)' freq(2*n+1:3*n)'];
freq_mean_std = [mean(freq_pre_tar_pos,1);std(freq_pre_tar_pos,0,1)];
%[124.530	62.216	3619.785
% 310.919	76.996	6725.196]

len_pre_tar_pos = [len(1:n)' len(n+1:2*n)' len(2*n+1:3*n)'];
len_mean_std = [mean(len_pre_tar_pos,1);std(len_pre_tar_pos,0,1)];
%[5.538 5.281	5.406
% 1.063	1.083	2.007]


%% location of the target words
tar_loc_mean_std = [mean(target_loc,1) std(target_loc,0,1)]; 
% [6.43 1.30]


%% number of words in the sentences
wrdnunm = cellfun(@isempty,sv_mat,'Uni',true);
wrdnunm = sum(~wrdnunm,2);

wrdnunm_mean_std = [mean(wrdnunm,1) std(wrdnunm,0,1)]; 
% [11.6 1.7]


%%
min(len_pre_tar_pos) %[4 4 2]
max(len_pre_tar_pos) %[8 7 12]

preAndtar = len_pre_tar_pos(:,1)+len_pre_tar_pos(:,2);
[min(preAndtar) max(preAndtar)] % [8 15]
[mean(preAndtar,1) std(preAndtar,0,1)] % [10.8 1.5]









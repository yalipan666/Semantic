%% 20220203 
% try to find correlations between tagging response, N400 amplitude
% difference, and eye movements


%% paths
Path.Tag = 'Z:\Semantic\Results\sv\Coh\';
Path.N400 = 'Z:\Semantic\Results\sv\ERF_allfixations\';
Path.EM = 'Z:\Semantic\Results\sv\Behav\';
Path.Result = 'Z:\Semantic\Results\sv\Corr\';
Corr = [];
Corr.CondName = {'Incong','Cong'};
    
%% get the averaged Tagging coherence for each participant
load([Path.Tag 'TagCoh']);    
Corr.TagSigSub = TagCoh.SigSubID;
Corr.TagPre = TagCoh.PreTarg_Ttest.data4test;


%% get the averaged N400 from sig sensors over sig time windows for each participant
load([Path.N400 'ERF_rmbl_cmb']);    
EpochType = {'PreTarg';'Targ';'PosTarg'};
n_sub = length(ERF.subs);
erf_tim = ERF.PreTarg_incong{1,1}.time';
for eee = 2:3
    eval(['stat = ERF.' EpochType{eee} '_stat;']);
    sig_tp = sum(stat.mask,1)~=0;
    SigTim = stat.time(sig_tp); % unit in s
    SigTim = dsearchn(erf_tim,SigTim(1)):dsearchn(erf_tim,SigTim(end));
    sigchanid = sum(stat.negclusterslabelmat(:,sig_tp),2)~=0;
    SigSen = stat.label(sigchanid);
    erf_avg = nan(n_sub,2);
    for s = 1:n_sub
        eval(['label_all = ERF.' EpochType{eee} '_incong{1,s}.label;']);
        sigchan_id = cellfun(@(x) find(strcmp(x,label_all)),SigSen,'Uni',true);
        eval(['tmp_1 = mean(mean(ERF.' EpochType{eee} '_incong{1,s}.avg(sigchan_id,SigTim),1),2);']);
        eval(['tmp_2 = mean(mean(ERF.' EpochType{eee} '_cong{1,s}.avg(sigchan_id,SigTim),1),2);']);
        erf_all(s,:) = [tmp_1 tmp_2];
    end
    eval(['Corr.N400' EpochType{eee} '=erf_all.*10e12;']);
end
clear ERF


%% get the eye movements metrics for each participant
load([Path.EM 'BehaData']);    
Corr.FirstFix_Targ = BehaData.FirstFix.Tag;
Corr.Regres_Targ = BehaData.Regression_percent.Tag;
% reading time for the whole sentence
DurationPerwrd = nan(length(TagCoh.subs),1); %averaged fixation duration of all words
Duration_wrdBefoPre = nan(length(TagCoh.subs),1); % averaged fixation duration of the beginning words that before pre-target
Duration_wrdFromPre = nan(length(TagCoh.subs),1); % averaged fixation duration of the beginning words that from pre-target
for s = 1:length(TagCoh.subs)
    load(['Z:\Semantic\Analyse_data\sv_' TagCoh.subs{s} '\EyeData.mat']);
    load(['Z:\Semantic\RawData\PTB_data\' TagCoh.subs{s} '.mat'],'Para');
    % id for sv_sentences
    sv_id = ismember(Para.CondMat,[11 12]); 
    
    %%% averaged fixation durations of all words 
    tmp_allfixdu = cellfun(@(x) sum(x(~isnan(x(:,7)),3)),EyeData.TrlFixdata,'Uni',true); %%[allfix_duration]
    % word number in all sentences
    wrdnum = sum(~cellfun(@isempty,Para.SentMat,'Uni',true),2);
    % divided by the number of words in each sentence 
    tmp = tmp_allfixdu./wrdnum;
    DurationPerwrd(s,1) = mean(tmp(sv_id));
    
    %%% averaged fixation duration of the beginning words that before pre-target
    targid = arrayfun(@(x) {x},Para.FlkWords(:,1),'Uni',true);
    tmp_allfixdu = cellfun(@(x,y) sum(x(x(:,7)<=y-2,3))/(y-2),EyeData.TrlFixdata,targid,'Uni',true); 
    Duration_wrdBefoPre(s,1) = mean(tmp_allfixdu(sv_id));
    
    %%% averaged fixation duration of the beginning words that from pre-target
    tmp_allfixdu = cellfun(@(x,y) sum(x(x(:,7)>=y-1,3)),EyeData.TrlFixdata,targid,'Uni',true); 
    tmp = tmp_allfixdu./(wrdnum-Para.FlkWords(:,1)+2);
    Duration_wrdFromPre(s,1) = mean(tmp(sv_id));
end
BehaData.DurationPerwrd = DurationPerwrd;
save([Path.EM 'BehaData'],'BehaData');
Corr.DurationPerwrd = DurationPerwrd;

  
%% Corr of Tagging and N400
% Tag_pretarg & N400_targ
tag_dif = Corr.TagPre(:,1)-Corr.TagPre(:,2);
N400_targ = Corr.N400Targ;
n400_dif = N400_targ(:,1)-N400_targ(:,2);
[coef,pval] = corr([tag_dif n400_dif(Corr.TagSigSub)],'type','Pearson') % p = 0.485
Corr.TagPre_N400Targ_CoefP = [coef(1,2) pval(1,2)];

% Tag_pretarg & N400_postarg
N400_postarg = Corr.N400_PosTarg(Corr.TagSigSub,:);
n400_dif = N400_postarg(:,1)-N400_postarg(:,2);
[coef, pval] = corr([tag_dif n400_dif],'type','Pearson')% p = 0.849


%% Corr of Tagging and EM
% Tag_pretarg & FirstFix_targ
ff = Corr.FirstFix_Targ;
ff_dif = ff(:,1)-ff(:,2);
[coef,pval] = corr([tag_dif ff_dif(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_FFdifTarg_CoefP = [coef(1,2),pval(1,2)]; % p = 0.631
ff_avg = mean(ff,2);
[coef,pval] = corr([tag_dif ff_avg(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_FFavgTarg_CoefP = [coef(1,2),pval(1,2)]; % p = 0.077
% Tag_pretarg & Gaze_targ
Corr.Gaze_Targ = BehaData.Gaze.Tag;
gz_avg = mean(Corr.Gaze_Targ,2);
[coef,pval] = corr([tag_dif gz_avg(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_GZavgTarg_CoefP = [coef(1,2),pval(1,2)]; % p = 0.0396

% Tag_pretarg & reading speed
rs = 1000./Corr.DurationPerwrd;
[coef,pval] = corr([tag_dif rs(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_ReadSped_CoefP = [coef(1,2),pval(1,2)]; % p = 0.048
[coef,pval] = corr([tag_dif Corr.DurationPerwrd(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_DurationPerwrd_CoefP = [coef(1,2),pval(1,2)]; % p = 0.0049

% Tag_pretarg & Reg_targ
reg = Corr.Regres_Targ;
reg_dif = reg(:,1)-reg(:,2);
[coef,pval] = corr([tag_dif reg_dif(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_RegdifTarg_CoefP = [coef(1,2),pval(1,2)]; % p = 0.714
reg_avg = mean(reg,2);
[coef,pval] = corr([tag_dif reg_avg(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_RegavgTarg_CoefP = [coef(1,2),pval(1,2)]; % p = 0.234

% Tag_pretarg & FirstFix_pretarg
[coef,pval] = corr([tag_dif mean(BehaData.FirstFix.Pre(Corr.TagSigSub,:),2)],'type','Pearson'); 
Corr.TagPre_FFavgPreTarg_CoefP = [coef(1,2),pval(1,2)]; %p = 0.024

% Tag_pretarg & FirstFix_targ_postarg
ff_t_p = mean(BehaData.FirstFix.Tag + BehaData.FirstFix.Pos,2);
[coef,pval] = corr([tag_dif ff_t_p(Corr.TagSigSub,:)],'type','Pearson'); 
Corr.TagPre_FFavgTargPostarg_CoefP = [coef(1,2),pval(1,2)]; %p = 0.0179


%% Corr of N400_targ and EM
% N400_targ & reading speed
[coef,pval] = corr([n400_dif Corr.DurationPerwrd],'type','Pearson');
Corr.N400Targ_DurationPerwrd_CoefP = [coef(1,2),pval(1,2)]; % p = 0.492

% N400_targ & Reg_targ
[coef,pval] = corr([n400_dif reg_dif],'type','Pearson')
Corr.N400Targ_RegdifTarg_CoefP = [coef(1,2),pval(1,2)]; % p = 0.0021

% N400_targ & FirstFix_targ_postarg
[coef,pval] = corr([n400_dif ff_t_p],'type','Pearson'); 
Corr.N400Targ_FFavgTargPostarg_CoefP = [coef(1,2),pval(1,2)]; %p = 0.1367

% N400_targ & Gaze_targ
[coef,pval] = corr([n400_dif mean(Corr.Gaze_Targ,2)],'type','Pearson'); 
Corr.N400Targ_GZavgTarg_CoefP = [coef(1,2),pval(1,2)]; %p = 0.240



%% Corr of EM: DurationPerwrd reg_dif
[coef,pval] = corr([DurationPerwrd reg_dif],'type','Pearson') % p = 0.0281
Corr.DurationPerwrd_RegdifTarg_CoefP = [coef(1,2),pval(1,2)]; 
save([Path.Result 'Corr'],'Corr')

%% simple plot
sigsub = Corr.TagSigSub;
Brain = {'TagPre','N400Targ'};
EM = {'DurationPerwrd','RegdifTarg'};
RegdifTarg = Corr.Regres_Targ(:,1)-Corr.Regres_Targ(:,2);
Corr.RegdifTarg = RegdifTarg;

for bb = 1:length(Brain)
    for em = 1:length(EM)
        eval(['b_tmp = Corr.' Brain{bb} ';']);
        b_tmp = b_tmp(:,1)-b_tmp(:,2);
        eval(['EM_tmp = Corr.' EM{em} ';']);
        if bb == 1
           EM_tmp = EM_tmp(sigsub);
        end
        eval(['cp = Corr.' Brain{bb} '_' EM{em} '_CoefP;']);
        cp = round(cp*1000)./1000;
        figtitle = [Brain{bb} '-' EM{em} ', c=' num2str(cp(1)) ', p=' num2str(cp(2))];
        h = figure('color',[1 1 1]);
        scatter(b_tmp,EM_tmp);
        xlabel(Brain{bb},'FontSize',12)
        ylabel(EM{em},'FontSize',12)
        title(figtitle,'FontSize',14)
        saveas(h,[Path.Result Brain{bb} '_' EM{em}]);
    end
end

%
b_tmp = Corr.TagPre(:,1)-Corr.TagPre(:,2);
EM_tmp = Corr.N400Targ(:,1)-Corr.N400Targ(:,2);
cp = Corr.TagPre_N400Targ_CoefP;
cp = round(cp*1000)./1000;
figtitle = ['TagPre-N400Targ, c=' num2str(cp(1)) ', p=' num2str(cp(2))];
h = figure('color',[1 1 1]);
scatter(b_tmp,EM_tmp(sigsub));
xlabel('TagPre','FontSize',12)
ylabel('N400Targ','FontSize',12)
title(figtitle,'FontSize',14)
saveas(h,[Path.Result 'TagPre_N400Targ']);






























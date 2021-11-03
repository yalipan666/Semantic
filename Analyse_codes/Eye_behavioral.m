%%% analyse of behavioral eye movement data
%%% edit: YPan, 20210918


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%========== stat on behavioral data ========== %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% gaze duration: all fixation on target before eye moves out of taget,
%re-read is not included
clear; clc;
load(['Z:\Analyse_data\SubInfo.mat']);
subjects = SubInfo.subjects.Targ60;
BehaData = [];
for sss = 1:length(subjects)
    load(['Z:\Analyse_data\' subjects{sss} '\Event.mat']);
    tmp = Event.event_raw; %%{'sentence_id'},{'word_loc'},{'loc2targ'},{'word_freq'},{'word_length'},{'saccade2this_du…'},{'fixation_on_MEG'},{'fixation_duration'},{'NextOrder'},{'FirstPassFix'},{'PreviousOrder'},{'SentenceCond'},{'PupilSize'}
    %%% remove outliers
    tmp = tmp(tmp(:,8)>=80 & tmp(:,8)<=1000,:);
    %%% get duration & pupil size: [low high]
    pre1 = tmp(tmp(:,3)==-1 & tmp(:,10)==1 & tmp(:,12)==1,[8 13]);
    pre2 = tmp(tmp(:,3)==-1 & tmp(:,10)==1 & tmp(:,12)==2,[8 13]);
    tag1 = tmp(tmp(:,3)== 0 & tmp(:,10)==1 & tmp(:,12)==1,[8 13]);
    tag2 = tmp(tmp(:,3)== 0 & tmp(:,10)==1 & tmp(:,12)==2,[8 13]);
    pos1 = tmp(tmp(:,3)== 1 & tmp(:,10)==1 & tmp(:,12)==1,[8 13]);
    pos2 = tmp(tmp(:,3)== 1 & tmp(:,10)==1 & tmp(:,12)==2,[8 13]);
    %%% fixation duration
    FirstFix.Pre(sss,:) = [mean(pre1(:,1))   mean(pre2(:,1))];
    FirstFix.Tag(sss,:) = [mean(tag1(:,1))   mean(tag2(:,1))];
    FirstFix.Pos(sss,:) = [mean(pos1(:,1))   mean(pos2(:,1))];
    %%% gaze and total-gaze, the saccades duration is not included here
    Gaze.Pre(sss,:) = [sum(tmp(tmp(:,3)==-1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==1,8))./size(pre1,1)  sum(tmp(tmp(:,3)==-1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==2,8))./size(pre2,1)];
    Gaze.Tag(sss,:) = [sum(tmp(tmp(:,3)== 0 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==1,8))./size(tag1,1)  sum(tmp(tmp(:,3)== 0 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==2,8))./size(tag2,1)];
    Gaze.Pos(sss,:) = [sum(tmp(tmp(:,3)== 1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==1,8))./size(pos1,1)  sum(tmp(tmp(:,3)== 1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==2,8))./size(pos2,1)];
    TotalGaze.Pre(sss,:) = [sum(tmp(tmp(:,3)==-1 & tmp(:,12)==1,8))./size(pre1,1)   sum(tmp(tmp(:,3)==-1 & tmp(:,12)==2,8))./size(pre2,1)];
    TotalGaze.Tag(sss,:) = [sum(tmp(tmp(:,3)== 0 & tmp(:,12)==1,8))./size(tag1,1)   sum(tmp(tmp(:,3)== 0 & tmp(:,12)==2,8))./size(tag2,1)];
    TotalGaze.Pos(sss,:) = [sum(tmp(tmp(:,3)== 1 & tmp(:,12)==1,8))./size(pos1,1)   sum(tmp(tmp(:,3)== 1 & tmp(:,12)==2,8))./size(pos2,1)];
    %%% pupil size of the first fixation
    PupilSize.Pre(sss,:) = [mean(pre1(:,2))   mean(pre2(:,2))];
    PupilSize.Tag(sss,:) = [mean(tag1(:,2))   mean(tag2(:,2))];
    PupilSize.Pos(sss,:) = [mean(pos1(:,2))   mean(pos2(:,2))];
end
[~,p] = ttest(FirstFix.Pre(:,1),FirstFix.Pre(:,2))
FirstFix.Pre_pvalue = p;
[~,p] = ttest(FirstFix.Tag(:,1),FirstFix.Tag(:,2))
FirstFix.Tag_pvalue = p;
[~,p] = ttest(FirstFix.Pos(:,1),FirstFix.Pos(:,2))
FirstFix.Pos_pvalue = p;

[~,p] = ttest(Gaze.Pre(:,1),Gaze.Pre(:,2))
Gaze.Pre_pvalue = p;
[~,p] = ttest(Gaze.Tag(:,1),Gaze.Tag(:,2))
Gaze.Tag_pvalue = p;
[~,p] = ttest(Gaze.Pos(:,1),Gaze.Pos(:,2))
Gaze.Pos_pvalue = p;

[~,p] = ttest(TotalGaze.Pre(:,1),TotalGaze.Pre(:,2))
TotalGaze.Pre_pvalue = p;
[~,p] = ttest(TotalGaze.Tag(:,1),TotalGaze.Tag(:,2))
TotalGaze.Tag_pvalue = p;
[~,p] = ttest(TotalGaze.Pos(:,1),TotalGaze.Pos(:,2))
TotalGaze.Pos_pvalue = p;

[~,p] = ttest(PupilSize.Pre(:,1),PupilSize.Pre(:,2))
PupilSize.Pre_pvalue = p;
[~,p] = ttest(PupilSize.Tag(:,1),PupilSize.Tag(:,2))
PupilSize.Tag_pvalue = p;
[~,p] = ttest(PupilSize.Pos(:,1),PupilSize.Pos(:,2))
PupilSize.Pos_pvalue = p;

BehaData.FirstFix = FirstFix;
BehaData.Gaze = Gaze;
BehaData.TotalGaze = TotalGaze;
BehaData.PupilSize = PupilSize;
save('Z:\Results\Behavioral\Targ60_BehaData.mat','BehaData');


%% %%%%%%%%%%%%%=============== plotting =================%%%%%%%%%%%%%%%%
FigNam = {'FirstFix';'Gaze';'TotalGaze';'PupilSize'};
for ff = 1:length(FigNam)
    eval(['tmp = BehaData.' FigNam{ff} ';']);
    h = figure('Name',FigNam{ff},'color',[1 1 1]);
    y = [mean(tmp.Pre); mean(tmp.Tag); mean(tmp.Pos)];         % random y values (3 groups of 2 parameters)
    errY = [std(tmp.Pre); std(tmp.Tag); std(tmp.Pos)]./sqrt(size(tmp.Pre,1));
    a = barwitherr(cat(3,-errY,errY), y);    % Plot with errorbars
    if ff == 1
        ylim([160 270])
        yyy = 245;
    elseif ff == 2
        ylim([200 350])
        yyy = 320;
    elseif ff == 3
        ylim([200 440])
        yyy = 395;
    else
        ylim([3500 4300])
        yyy = 4150;
    end
    set(gca,'box','off','LineWidth',2)
    set(gca,'XTickLabel',{'Pre-target','Target','Post-target'},'FontSize',16,'FontWeight','bold')
    text(1,yyy,num2str(round(1000*tmp.Pre_pvalue)/1000),'FontSize',12,'FontWeight','bold');
    text(2,yyy,num2str(round(1000*tmp.Tag_pvalue)/1000),'FontSize',12,'FontWeight','bold');
    text(3,yyy,num2str(round(1000*tmp.Pos_pvalue)/1000),'FontSize',12,'FontWeight','bold');
    switch ff
        case 1
            subtitle = 'First Fixation duration(ms)';
        case 2
            subtitle = 'Gaze duration(ms)';
        case 3
            subtitle = 'Total Gaze duration(ms)';
        case 4
            subtitle = 'Pupil Dilation';
    end
    title(subtitle,'FontSize',24,'FontWeight','bold');
    legendflex(a,[{'low'},{'high'}],'anchor', {'ne','ne'}, 'buffer', [5 -5],'fontsize',12,'FontWeight','bold','xscale',0.8,'box', 'off');
    saveas(h,['Z:\Results\Behavioral\Targ60_' FigNam{ff}]);
end



% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%==========JEP ========== %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load('U:\ReadingProject\stimuli\JEP_stimuli\WordsFrequencyStat_JEP.mat'); % AllWords
%%% version 1 (exact the same for version 2, but reversed for conditions)
load('U:\ReadingProject\PTB_codes\Reading_Exp\TarWodLocFreq_JEP_1.mat');
TargLoc = TarWodLocFreq(:,[1 3]);
Targ = []; preTarg = []; posTarg = [];
for tt = 1:length(TargLoc)
    Targ(tt,:)    = [AllWords.V1_celex(tt,TargLoc(tt,1)) AllWords.V1_length(tt,TargLoc(tt,1)) AllWords.V1_celex(tt,TargLoc(tt,2)) AllWords.V1_length(tt,TargLoc(tt,2))]; %[WordFreq WordLength]
    preTarg(tt,:) = [AllWords.V1_celex(tt,TargLoc(tt,1)-1) AllWords.V1_length(tt,TargLoc(tt,1)-1) AllWords.V1_celex(tt,TargLoc(tt,2)-1) AllWords.V1_length(tt,TargLoc(tt,2)-1)]; %[WordFreq WordLength]
    posTarg(tt,:) = [AllWords.V1_celex(tt,TargLoc(tt,1)+1) AllWords.V1_length(tt,TargLoc(tt,1)+1) AllWords.V1_celex(tt,TargLoc(tt,2)+1) AllWords.V1_length(tt,TargLoc(tt,2)+1)]; %[WordFreq WordLength]
end
Targ = [Targ(:,[1 2]); Targ(:,[3 4]);];
preTarg = [preTarg(:,[1 2]); preTarg(:,[3 4]);];
posTarg = [posTarg(:,[1 2]); posTarg(:,[3 4]);];
freqid = [TarWodLocFreq(:,2);TarWodLocFreq(:,2)];
idx_targ_low = find(freqid == 1);
idx_targ_high = find(freqid == 2);
Targ_low = Targ(idx_targ_low,:);
Targ_high = Targ(idx_targ_high,:);
preTarg_low = preTarg(idx_targ_low,:);
preTarg_high = preTarg(idx_targ_high,:);
posTarg_low = posTarg(idx_targ_low,:);
posTarg_high = posTarg(idx_targ_high,:);
%%% [target_word_frequency  word_length] effect on target words
%%% ttest for both frequency and length--p = [freq length]
[~,p1] = ttest(Targ_low,Targ_high)  % p1 = 0.0000    0.0416
[~,p2] = ttest(preTarg_low,preTarg_high) % p2 = 0.3215    0.8812
[~,p3] = ttest(posTarg_low,posTarg_high) % p3 = 0.0005    0.0056
%%% value:[freq-low length-low freq-high length-high]
nanmean([Targ_low,Targ_high])       %%%[5.3500    5.6744  109.7784  5.5233]
nanmean([preTarg_low,preTarg_high]) %%%[731.3788  5.3372  514.6147  5.3605]
nanmean([posTarg_low,posTarg_high]) %%%[1611.6    5.5233  402.16    6.1744]

%%% conlusion: Targets: sig effect on frequency but also length
%%% no sig effect on neith freq nor length for Pretargs, but both are sig
%%% for post-targets
%%% low freq targets have longer length

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%========== 2. stat on behavioral data ========== %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% gaze duration: all fixation on target before eye move out of taget,
%%% re-read is not included
clear
load(['Z:\Analyse_data\SubInfo.mat']);
subjects = SubInfo.subjects.JEP60;
% subjects = subjects(1:24); %%% equality version 1 and version 2
%%%% combining the two target data
BehaData = [];
for sss = 1:length(subjects)
    load(['Z:\Analyse_data\' subjects{sss} '\Event.mat']);
    tmp = Event.event_raw; %%{'sentence_id'},{'word_loc'},{'loc2targ'},{'word_freq'},{'word_length'},{'saccade2this_du…'},{'fixation_on_MEG'},{'fixation_duration'},{'NextOrder'},{'FirstPassFix'},{'PreviousOrder'},{SentenceCond}
    %%% remove outliers
    tmp = tmp(tmp(:,8)>=80 & tmp(:,8)<=1000 ,:);
    %%% get duration & pupil size: [low high]
    pre1 = tmp(tmp(:,3)==-1 & tmp(:,10)==1 & tmp(:,12)==1,[8 13]);
    pre2 = tmp(tmp(:,3)==-1 & tmp(:,10)==1 & tmp(:,12)==2,[8 13]);
    tag1 = tmp(tmp(:,3)== 0 & tmp(:,10)==1 & tmp(:,12)==1,[8 13]);
    tag2 = tmp(tmp(:,3)== 0 & tmp(:,10)==1 & tmp(:,12)==2,[8 13]);
    pos1 = tmp(tmp(:,3)== 1 & tmp(:,10)==1 & tmp(:,12)==1,[8 13]);
    pos2 = tmp(tmp(:,3)== 1 & tmp(:,10)==1 & tmp(:,12)==2,[8 13]);
    %%% fixation duration
    FirstFix.Pre(sss,:) = [mean(pre1(:,1))   mean(pre2(:,1))];
    FirstFix.Tag(sss,:) = [mean(tag1(:,1))   mean(tag2(:,1))];
    FirstFix.Pos(sss,:) = [mean(pos1(:,1))   mean(pos2(:,1))];
    Gaze.Pre(sss,:) = [sum(tmp(tmp(:,3)==-1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==1,8))./size(pre1,1)   sum(tmp(tmp(:,3)==-1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==2,8))./size(pre2,1)];
    Gaze.Tag(sss,:) = [sum(tmp(tmp(:,3)== 0 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==1,8))./size(tag1,1)  sum(tmp(tmp(:,3)== 0 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==2,8))./size(tag2,1)];
    Gaze.Pos(sss,:) = [sum(tmp(tmp(:,3)== 1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==1,8))./size(pos1,1)   sum(tmp(tmp(:,3)== 1 & (tmp(:,10)==1 | tmp(:,11)==0) & tmp(:,12)==2,8))./size(pos2,1)];
    TotalGaze.Pre(sss,:) = [sum(tmp(tmp(:,3)==-1 & tmp(:,12)==1,8))./size(pre1,1)   sum(tmp(tmp(:,3)==-1 & tmp(:,12)==2,8))./size(pre2,1)];
    TotalGaze.Tag(sss,:) = [sum(tmp(tmp(:,3)== 0 & tmp(:,12)==1,8))./size(tag1,1)   sum(tmp(tmp(:,3)== 0 & tmp(:,12)==2,8))./size(tag2,1)];
    TotalGaze.Pos(sss,:) = [sum(tmp(tmp(:,3)== 1 & tmp(:,12)==1,8))./size(pos1,1)   sum(tmp(tmp(:,3)== 1 & tmp(:,12)==2,8))./size(pos2,1)];
    %%% pupil size
    PupilSize.Pre(sss,:) = [mean(pre1(:,2))   mean(pre2(:,2))];
    PupilSize.Tag(sss,:) = [mean(tag1(:,2))   mean(tag2(:,2))];
    PupilSize.Pos(sss,:) = [mean(pos1(:,2))   mean(pos2(:,2))];
end
[~,p] = ttest(FirstFix.Pre(:,1),FirstFix.Pre(:,2))
FirstFix.Pre_pvalue = p;
[~,p] = ttest(FirstFix.Tag(:,1),FirstFix.Tag(:,2))
FirstFix.Tag_pvalue = p;
[~,p] = ttest(FirstFix.Pos(:,1),FirstFix.Pos(:,2))
FirstFix.Pos_pvalue = p;

[~,p] = ttest(Gaze.Pre(:,1),Gaze.Pre(:,2))
Gaze.Pre_pvalue = p;
[~,p] = ttest(Gaze.Tag(:,1),Gaze.Tag(:,2))
Gaze.Tag_pvalue = p;
[~,p] = ttest(Gaze.Pos(:,1),Gaze.Pos(:,2))
Gaze.Pos_pvalue = p;

[~,p] = ttest(TotalGaze.Pre(:,1),TotalGaze.Pre(:,2))
TotalGaze.Pre_pvalue = p;
[~,p] = ttest(TotalGaze.Tag(:,1),TotalGaze.Tag(:,2))
TotalGaze.Tag_pvalue = p;
[~,p] = ttest(TotalGaze.Pos(:,1),TotalGaze.Pos(:,2))
TotalGaze.Pos_pvalue = p;

[~,p] = ttest(PupilSize.Pre(:,1),PupilSize.Pre(:,2));
PupilSize.Pre_pvalue = p;
[~,p] = ttest(PupilSize.Tag(:,1),PupilSize.Tag(:,2))
PupilSize.Tag_pvalue = p;
[~,p] = ttest(PupilSize.Pos(:,1),PupilSize.Pos(:,2))
PupilSize.Pos_pvalue = p;


BehaData.FirstFix = FirstFix;
BehaData.Gaze = Gaze;
BehaData.TotalGaze = TotalGaze;
BehaData.PupilSize = PupilSize;

save('Z:\Results\Behavioral\JEP60_BehaData.mat','BehaData');

%% %%%%%%%%%%%%%=============== plotting =================%%%%%%%%%%%%%%%%
FigNam = {'FirstFix';'Gaze';'TotalGaze';}; %'PupilSize'};
figtitle = 'JEP_BehaData';
h = figure('Name',figtitle,'color',[1 1 1]);
ScSz = [1 1 1000 400];%get(groot, 'Screensize' );
set(h,'Position', ScSz,'color',[1 1 1],'MenuBar','figure');
pos_grid =  [ScSz(1)+20 ScSz(2)+20  ScSz(3)-20  ScSz(4)-20];
nrow = 1;
ncol = length(FigNam);
gridPos = DivideScreen(nrow,ncol,ScSz,90,110,pos_grid);
for ff = 1:length(FigNam)
    loc = ff;
    h = axes('Position',gridPos{loc});
    eval(['tmp = BehaData.' FigNam{ff} ';']);
    y = [mean(tmp.Pre); mean(tmp.Tag); mean(tmp.Pos)];         % random y values (3 groups of 2 parameters)
    errY = [std(tmp.Pre); std(tmp.Tag); std(tmp.Pos)]./sqrt(size(tmp.Pre,1));
    a = barwitherr(cat(3,-errY,errY), y);    % Plot with errorbars
    if ff == 1
        ylim([170 250])
        yyy = 250;
    elseif ff == 2
        ylim([200 300])
        yyy = 280;
    elseif ff == 3
        ylim([250 350])
        yyy = 350;
    else
        ylim([3500 4100])
        yyy = 4000;
    end
    set(gca,'box','off','LineWidth',2)
    set(gca,'XTickLabel',{'PreTarget','Target','PosTarget'},'FontSize',12,'FontWeight','bold')
    text(1,yyy,num2str(round(1000*tmp.Pre_pvalue)/1000),'FontSize',12,'FontWeight','bold');
    text(2,yyy,num2str(round(1000*tmp.Tag_pvalue)/1000),'FontSize',12,'FontWeight','bold');
    text(3,yyy,num2str(round(1000*tmp.Pos_pvalue)/1000),'FontSize',12,'FontWeight','bold');
    switch ff
        case 1
            subtitle = 'First Fixation (ms)';
        case 2
            subtitle = 'Gaze (ms)';
        case 3
            subtitle = 'Total Gaze (ms)';
        case 4
            subtitle = 'Pupil Dilation';
    end
    title(subtitle,'FontSize',16,'FontWeight','bold');
    legendflex(a,[{'low'},{'high'}],'anchor', {'ne','ne'}, 'buffer', [5 -5],'fontsize',12,'FontWeight','bold','xscale',0.8,'box', 'off');
    axcopy;
end
saveas(h,['Z:\Results\Behavioral\' figtitle]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%========= full statistics on word properties=========%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get statistics on word length
%%% divided data acorrding to the word length
if exist('Z:\Analyse_data\WordLength.mat','file')
    load('Z:\Analyse_data\WordLength.mat')
else
    WordLength = [];
end
DataSets = {'Targ60';'Targ5560';'PosTarg60';'JEP60'};

%% %%=============Targ60 with 2 setence versions
dd = 1;
sentversion = '_2';
load('U:\ReadingProject\stimuli\WordsFrequencyStat.mat'); % AllWords
load(['U:\ReadingProject\PTB_codes\Reading_Exp\TarWodLocFreq' sentversion '.mat']); %%V1, same word length for pretarget-target-posttarget for both versions in JEP, but not Targ60
eval(['WrdLen_all = AllWords.V' sentversion(2) '_length;']);
eval(['WrdLen_freq = AllWords.V' sentversion(2) '_celex;']);
n = size(WrdLen_all,1);
WrdLenFreq_PreTarg = [zeros(n,1) TarWodLocFreq(:,2) zeros(n,1)];
WrdLenFreq_Targ = [zeros(n,1) TarWodLocFreq(:,2) zeros(n,1)];
WrdLenFreq_PosTarg = [zeros(n,1) TarWodLocFreq(:,2) zeros(n,1)];
for ww = 1:n
    WrdLenFreq_PreTarg(ww,[1 3]) = [WrdLen_all(ww,TarWodLocFreq(ww,1)-1) WrdLen_freq(ww,TarWodLocFreq(ww,1)-1)];
    WrdLenFreq_Targ(ww,[1 3])    = [WrdLen_all(ww,TarWodLocFreq(ww,1))   WrdLen_freq(ww,TarWodLocFreq(ww,1))];
    WrdLenFreq_PosTarg(ww,[1 3]) = [WrdLen_all(ww,TarWodLocFreq(ww,1)+1) WrdLen_freq(ww,TarWodLocFreq(ww,1)+1)];
end
%%% get the word length that in the middle of distribution
WrdLen_PreTarg_tmp = sort(WrdLenFreq_PreTarg(:,1)); mid_PreTarg = WrdLen_PreTarg_tmp(n/2);
WrdLen_Targ_tmp = sort(WrdLenFreq_Targ(:,1)); mid_Targ = WrdLen_Targ_tmp(n/2);
WrdLen_PosTarg_tmp = sort(WrdLenFreq_PosTarg(:,1)); mid_PosTarg = WrdLen_PosTarg_tmp(n/2);
eval(['WordLength.' DataSets{dd} sentversion '.WrdLen_PreTarg = WrdLenFreq_PreTarg;'])
eval(['WordLength.' DataSets{dd} sentversion '.WrdLen_Targ = WrdLenFreq_Targ;'])
eval(['WordLength.' DataSets{dd} sentversion '.WrdLen_PosTarg = WrdLenFreq_PosTarg;'])
eval(['WordLength.' DataSets{dd} sentversion '.MidNum_PreTarg = mid_PreTarg;'])
eval(['WordLength.' DataSets{dd} sentversion '.MidNum_Targ = mid_Targ;'])
eval(['WordLength.' DataSets{dd} sentversion '.MidNum_PosTarg = mid_PosTarg;'])

%%%%%======= stat for Targ60: get conditions based on the target-length =======%%%%%
%%% get sentence id for conditions(long/short targets)
shrt_id = find(WrdLenFreq_Targ(:,1)<mid_Targ); %% 50 ---[5]
long_id = find(WrdLenFreq_Targ(:,1)>mid_Targ); %% 50 ---[7 8 9]
PreFreq_shrt = [];
PreFreq_long = [];
PreLength_shrt = [];
PreLength_long = [];
TargFreq_shrt = [];
TargFreq_long = [];
TargLength_shrt = [];
TargLength_long = [];
PosFreq_shrt = [];
PosFreq_long = [];
PosLength_shrt = [];
PosLength_long = [];
n = length(long_id); %%% length of long_id should be the same as short_id
for ww = 1:n
    PreFreq_shrt(ww,1) = WrdLen_freq(shrt_id(ww),TarWodLocFreq(shrt_id(ww),1)-1);
    PreFreq_long(ww,1) = WrdLen_freq(long_id(ww),TarWodLocFreq(long_id(ww),1)-1);
    PreLength_shrt(ww,1) = WrdLen_all(shrt_id(ww),TarWodLocFreq(shrt_id(ww),1)-1);
    PreLength_long(ww,1) = WrdLen_all(long_id(ww),TarWodLocFreq(long_id(ww),1)-1);
    TargFreq_shrt(ww,1) = WrdLen_freq(shrt_id(ww),TarWodLocFreq(shrt_id(ww),1));
    TargFreq_long(ww,1) = WrdLen_freq(long_id(ww),TarWodLocFreq(long_id(ww),1));
    TargLength_shrt(ww,1) = WrdLen_all(shrt_id(ww),TarWodLocFreq(shrt_id(ww),1));
    TargLength_long(ww,1) = WrdLen_all(long_id(ww),TarWodLocFreq(long_id(ww),1));
    PosFreq_shrt(ww,1) = WrdLen_freq(shrt_id(ww),TarWodLocFreq(shrt_id(ww),1)+1);
    PosFreq_long(ww,1) = WrdLen_freq(long_id(ww),TarWodLocFreq(long_id(ww),1)+1);
    PosLength_shrt(ww,1) = WrdLen_all(shrt_id(ww),TarWodLocFreq(shrt_id(ww),1)+1);
    PosLength_long(ww,1) = WrdLen_all(long_id(ww),TarWodLocFreq(long_id(ww),1)+1);
end
[~,p,~,stats] = ttest(PreFreq_shrt, PreFreq_long);
stats.p = p; disp(['PreTargFreq_CondTargLeng = ' num2str(p)])
stats.PreFreq_shrt_long = [PreFreq_shrt PreFreq_long];
eval(['WordLength.' DataSets{dd} sentversion '.Tstat.PreTargFreq_CondTargLeng = stats;'])
[~,p,~,stats] = ttest(PreLength_shrt, PreLength_long);
stats.p = p; disp(['PreTargLength_CondTargLeng = ' num2str(p)])
stats.PreLength_shrt_long = [PreLength_shrt PreLength_long];
eval(['WordLength.' DataSets{dd} sentversion '.Tstat.PreTargLength_CondTargLeng = stats;'])

[~,p,~,stats] = ttest(TargFreq_shrt, TargFreq_long);
stats.p = p; disp(['TargFreq_CondTargLeng = ' num2str(p)])
stats.TargFreq_shrt_long = [TargFreq_shrt TargFreq_long];
eval(['WordLength.' DataSets{dd} sentversion '.Tstat.TargFreq_CondTargLeng = stats;'])
[~,p,~,stats] = ttest(TargLength_shrt, TargLength_long);
stats.p = p; disp(['TargLength_CondTargLeng = ' num2str(p)])
stats.TargLength_shrt_long = [TargLength_shrt TargLength_long];
eval(['WordLength.' DataSets{dd} sentversion '.Tstat.TargLength_CondTargLeng = stats;'])

[~,p,~,stats] = ttest(PosFreq_shrt, PosFreq_long);
stats.p = p; disp(['PosTargFreq_CondTargLeng = ' num2str(p)])
stats.PosFreq_shrt_long = [PosFreq_shrt PosFreq_long];
eval(['WordLength.' DataSets{dd} sentversion '.Tstat.PosTargFreq_CondTargLeng = stats;'])
[~,p,~,stats] = ttest(PosLength_shrt, PosLength_long);
stats.p = p; disp(['PosTargLength_CondTargLeng = ' num2str(p)])
stats.PosLength_shrt_long = [PosLength_shrt PosLength_long];
eval(['WordLength.' DataSets{dd} sentversion '.Tstat.PosTargLength_CondTargLeng = stats;'])

mean([PreLength_long PreLength_shrt TargLength_long TargLength_shrt]) %%% [7.2800    6.7600    7.6400    5.0000]
save('Z:\Analyse_data\WordLength','WordLength');

%%% plot figure
h = figure('Name',['WordLength_' DataSets{dd}],'color',[1 1 1]);
hhh = subplot(3,3,1); hist(WrdLenFreq_PreTarg(:,1),20); title('PreTarg WordLength All');
hold on; plot([mid_PreTarg mid_PreTarg],get(hhh,'ylim'),'-.r','LineWidth',2)
subplot(3,3,4); hist(WrdLenFreq_PreTarg(WrdLenFreq_Targ(:,2)==1,1),20); title('PreTarg WordLength LowFreq');
subplot(3,3,7); hist(WrdLenFreq_PreTarg(WrdLenFreq_Targ(:,2)==2,1),20); title('PreTarg WordLength HighFreq');
hhh = subplot(3,3,2); hist(WrdLenFreq_Targ(:,1),20); title('Targ WordLength All');
hold on; plot([mid_Targ mid_Targ],get(hhh,'ylim'),'-.r','LineWidth',2)
subplot(3,3,5); hist(WrdLenFreq_Targ(WrdLenFreq_Targ(:,2)==1,1),20); title('Targ WordLength LowFreq');
subplot(3,3,8); hist(WrdLenFreq_Targ(WrdLenFreq_Targ(:,2)==2,1),20); title('Targ WordLength HighFreq');
hhh = subplot(3,3,3); hist(WrdLenFreq_PosTarg(:,1),20); title('PosTarg WordLength All');
hold on; plot([mid_PosTarg mid_PosTarg],get(hhh,'ylim'),'-.r','LineWidth',2)
subplot(3,3,6); hist(WrdLenFreq_PosTarg(WrdLenFreq_Targ(:,2)==1,1),20); title('PosTarg WordLength LowFreq');
subplot(3,3,9); hist(WrdLenFreq_PosTarg(WrdLenFreq_Targ(:,2)==2,1),20); title('PosTarg WordLength HighFreq');
saveas(h,['Z:\Results\Behavioral\WordLength_' DataSets{dd}]);


%% ==========================%%% for JEP sentences with 2 targets each sentence
dd = 4;
load('U:\ReadingProject\stimuli\JEP_stimuli\WordsFrequencyStat_JEP.mat'); % AllWords
load('U:\ReadingProject\PTB_codes\Reading_Exp\TarWodLocFreq_JEP_1.mat'); %% same for both versions
WrdLen_all = AllWords.V1_length;
WrdFreq_all = AllWords.V1_celex;
n = size(WrdLen_all,1)*2;
TarWodLocFreq = [TarWodLocFreq(:,[1 2]); TarWodLocFreq(:,[3 2])];
WrdLenFreq_PreTarg_2 = [zeros(n,1) TarWodLocFreq(:,2) zeros(n,1)];
WrdLenFreq_Targ_2 = [zeros(n,1) TarWodLocFreq(:,2) zeros(n,1)];
WrdLenFreq_PosTarg_2 = [zeros(n,1) TarWodLocFreq(:,2) zeros(n,1)];
for ww = 1:n
    sid = ww; %% sentence id
    if sid > size(WrdLen_all,1)
        sid = sid - size(WrdLen_all,1);
    end
    WrdLenFreq_PreTarg_2(ww,[1 3]) = [WrdLen_all(sid,TarWodLocFreq(ww,1)-1) WrdFreq_all(sid,TarWodLocFreq(ww,1)-1)];
    WrdLenFreq_Targ_2(ww,[1 3]) = [WrdLen_all(sid,TarWodLocFreq(ww,1)) WrdFreq_all(sid,TarWodLocFreq(ww,1))];
    WrdLenFreq_PosTarg_2(ww,[1 3]) = [WrdLen_all(sid,TarWodLocFreq(ww,1)+1) WrdFreq_all(sid,TarWodLocFreq(ww,1)+1)];
end
%%% get the word length that in the middle of distribution
WrdLen_PreTarg_tmp = sort(WrdLenFreq_PreTarg_2(:,1)); mid_PreTarg = WrdLen_PreTarg_tmp(n/2);
WrdLen_Targ_tmp = sort(WrdLenFreq_Targ_2(:,1)); mid_Targ = WrdLen_Targ_tmp(n/2);
WrdLen_PosTarg_tmp = sort(WrdLenFreq_PosTarg_2(:,1)); mid_PosTarg = WrdLen_PosTarg_tmp(n/2);
eval(['WordLength.' DataSets{dd} '.WrdLen_PreTarg = WrdLenFreq_PreTarg_2;'])
eval(['WordLength.' DataSets{dd} '.WrdLen_Targ = WrdLenFreq_Targ_2;'])
eval(['WordLength.' DataSets{dd} '.WrdLen_PosTarg = WrdLenFreq_PosTarg_2;'])
eval(['WordLength.' DataSets{dd} '.MidNum_PreTarg = mid_PreTarg;'])
eval(['WordLength.' DataSets{dd} '.MidNum_Targ = mid_Targ;'])
eval(['WordLength.' DataSets{dd} '.MidNum_PosTarg = mid_PosTarg;'])
save('Z:\Analyse_data\WordLength','WordLength');

%%% plot figure
h = figure('Name',['WordLength_' DataSets{dd}],'color',[1 1 1]);
hhh = subplot(3,3,1); hist(WrdLenFreq_PreTarg_2(:,1),20); title('PreTarg WordLength All');
hold on; plot([mid_PreTarg mid_PreTarg],get(hhh,'ylim'),'-.r','LineWidth',2)
subplot(3,3,4); hist(WrdLenFreq_PreTarg_2(WrdLenFreq_Targ_2(:,2)==1,1),20); title('PreTarg WordLength LowFreq');
subplot(3,3,7); hist(WrdLenFreq_PreTarg_2(WrdLenFreq_Targ_2(:,2)==2,1),20); title('PreTarg WordLength HighFreq');
hhh = subplot(3,3,2); hist(WrdLenFreq_Targ_2(:,1),20); title('Targ WordLength All');
hold on; plot([mid_Targ mid_Targ],get(hhh,'ylim'),'-.r','LineWidth',2)
subplot(3,3,5); hist(WrdLenFreq_Targ_2(WrdLenFreq_Targ_2(:,2)==1,1),20); title('Targ WordLength LowFreq');
subplot(3,3,8); hist(WrdLenFreq_Targ_2(WrdLenFreq_Targ_2(:,2)==2,1),20); title('Targ WordLength HighFreq');
hhh = subplot(3,3,3); hist(WrdLenFreq_PosTarg_2(:,1),20); title('PosTarg WordLength All');
hold on; plot([mid_PosTarg mid_PosTarg],get(hhh,'ylim'),'-.r','LineWidth',2)
subplot(3,3,6); hist(WrdLenFreq_PosTarg_2(WrdLenFreq_Targ_2(:,2)==1,1),20); title('PosTarg WordLength LowFreq');
subplot(3,3,9); hist(WrdLenFreq_PosTarg_2(WrdLenFreq_Targ_2(:,2)==2,1),20); title('PosTarg WordLength HighFreq');
saveas(h,['Z:\Results\Behavioral\WordLength_' DataSets{dd}]);


%%%%%======= stat for JEP: get conditions based on target-freq ========%%%%%
%%% word length difference for pre-target words between high/low target
%%% frequency
PreTarg_WrdFreq_LowCond = WrdLenFreq_PreTarg_2(WrdLenFreq_Targ_2(:,2)==1,3);
PreTarg_WrdFreq_HighCond = WrdLenFreq_PreTarg_2(WrdLenFreq_Targ_2(:,2)==2,3);
[~,p,~,stats] = ttest(PreTarg_WrdFreq_LowCond, PreTarg_WrdFreq_HighCond);
stats.p = p; disp(['PreTargFreq_CondTargFreq_T = ' num2str(p)])
eval(['WordLength.' DataSets{dd} '.PreTargFreq_CondTargFreq_T = stats;'])

PreTarg_WrdLeng_LowCond = WrdLenFreq_PreTarg_2(WrdLenFreq_Targ_2(:,2)==1,1);
PreTarg_WrdLeng_HighCond = WrdLenFreq_PreTarg_2(WrdLenFreq_Targ_2(:,2)==2,1);
[~,p,~,stats] = ttest(PreTarg_WrdLeng_LowCond, PreTarg_WrdLeng_HighCond);
stats.p = p; disp(['PreTargLength_CondTargFreq_T = ' num2str(p)])
eval(['WordLength.' DataSets{dd} '.PreTargLength_CondTargFreq_T = stats;'])

Targ_WrdFreq_LowCond = WrdLenFreq_Targ_2(WrdLenFreq_Targ_2(:,2)==1,3);
Targ_WrdFreq_HighCond = WrdLenFreq_Targ_2(WrdLenFreq_Targ_2(:,2)==2,3);
[~,p,~,stats] = ttest(Targ_WrdFreq_LowCond, Targ_WrdFreq_HighCond);
stats.p = p; disp(['TargFreq_CondTargFreq_T = ' num2str(p)])
eval(['WordLength.' DataSets{dd} '.TargFreq_CondTargFreq_T = stats;'])

Targ_WrdLeng_LowCond = WrdLenFreq_Targ_2(WrdLenFreq_Targ_2(:,2)==1,1);
Targ_WrdLeng_HighCond = WrdLenFreq_Targ_2(WrdLenFreq_Targ_2(:,2)==2,1);
[~,p,~,stats] = ttest(Targ_WrdLeng_LowCond, Targ_WrdLeng_HighCond);
stats.p = p; disp(['TargLength_CondTargLeng_T = ' num2str(p)])
eval(['WordLength.' DataSets{dd} '.TargLength_CondTargLeng_T = stats;'])

PosTarg_WrdFreq_LowCond = WrdLenFreq_PosTarg_2(WrdLenFreq_Targ_2(:,2)==1,3);
PosTarg_WrdFreq_HighCond = WrdLenFreq_PosTarg_2(WrdLenFreq_Targ_2(:,2)==2,3);
[~,p,~,stats] = ttest(PosTarg_WrdFreq_LowCond, PosTarg_WrdFreq_HighCond);
stats.p = p; disp(['PosTargFreq_CondTargFreq_T = ' num2str(p)])
eval(['WordLength.' DataSets{dd} '.PosTargFreq_CondTargFreq_T = stats;'])

PosTarg_WrdLeng_LowCond = WrdLenFreq_PosTarg_2(WrdLenFreq_Targ_2(:,2)==1,1);
PosTarg_WrdLeng_HighCond = WrdLenFreq_PosTarg_2(WrdLenFreq_Targ_2(:,2)==2,1);
[~,p,~,stats] = ttest(PosTarg_WrdLeng_LowCond, PosTarg_WrdLeng_HighCond);
stats.p = p; disp(['PosTargLength_CondTargLeng_T = ' num2str(p)])
eval(['WordLength.' DataSets{dd} '.PosTargLength_CondTargLeng_T = stats;'])
mean([Targ_WrdFreq_LowCond Targ_WrdLeng_LowCond Targ_WrdFreq_HighCond  Targ_WrdLeng_HighCond])
% % [5.3500    5.6744  109.7784    5.5233]
% % PreTargFreq_CondTargFreq_T = 0.32148
% % PreTargLength_CondTargFreq_T = 0.8812
% % TargFreq_CondTargFreq_T = 9.6337e-07
% % TargLength_CondTargLeng_T = 0.041622
% % PosTargFreq_CondTargFreq_T = 0.00046359
% % PosTargLength_CondTargLeng_T = 0.0055987
mean([PosTarg_WrdFreq_LowCond PosTarg_WrdFreq_HighCond PosTarg_WrdLeng_LowCond  PosTarg_WrdLeng_HighCond])
[1611.6    402.2    5.5    6.2]


%%%% combing two versions together to get the perfect controls
WrdLenFreq_PreTarg_2 = [WrdLenFreq_PreTarg_1; WrdLenFreq_PreTarg_2];
WrdLenFreq_Targ_2 = [WrdLenFreq_Targ_1; WrdLenFreq_Targ_2 ];
WrdLenFreq_PosTarg_2 = [WrdLenFreq_PosTarg_1; WrdLenFreq_PosTarg_2 ];
%%% only: TargFreq_CondTargFreq_T = 9.6824e-14, other ttest are all p=1


WrdLenFreq_PreTarg = WrdLenFreq_PreTarg_1;
WrdLenFreq_Targ = WrdLenFreq_Targ_1;
WrdLenFreq_PosTarg = WrdLenFreq_PosTarg_1;


%% %%%======= stat for JEP: get conditions based on target-length ========%%%%%
%%%% cannot do the same loop stat for targ60, since it has the same number
%%%% of trials in long-short-target-length-condition already!
%%% word length difference for pre-target words between high/low target
%%% length
%%% since the number of trials are different for long/short target (i.e.
%%% target-length condition!), using loops to get the averaged significance
loops = 10000;
P_all = zeros(loops,6);
tvalue_all = zeros(loops,6);
sd_all = zeros(loops,6);
longid_all = find(WrdLenFreq_Targ(:,1)==6);
shortid = find(WrdLenFreq_Targ(:,1)==5);
for i = 1:loops
    tmp_id = randperm(length(longid_all))';
    longid = tmp_id(1:length(shortid));
    
    PreTarg_WrdFreq_LongCond = WrdLenFreq_PreTarg(longid,3);
    PreTarg_WrdFreq_ShortCond = WrdLenFreq_PreTarg(shortid,3);
    [~,p,~,stats] = ttest(PreTarg_WrdFreq_LongCond, PreTarg_WrdFreq_ShortCond);
    P_all(i,1) = p; tvalue_all(i,1) = stats.tstat; sd_all(i,1) = stats.sd;
    
    PreTarg_WrdLeng_LongCond = WrdLenFreq_PreTarg(longid,1);
    PreTarg_WrdLeng_ShortCond = WrdLenFreq_PreTarg(shortid,1);
    [~,p,~,stats] = ttest(PreTarg_WrdLeng_LongCond, PreTarg_WrdLeng_ShortCond);
    P_all(i,2) = p; tvalue_all(i,2) = stats.tstat; sd_all(i,2) = stats.sd;
    
    Targ_WrdFreq_LongCond = WrdLenFreq_Targ(longid,3);
    Targ_WrdFreq_ShortCond = WrdLenFreq_Targ(shortid,3);
    [~,p,~,stats] = ttest(Targ_WrdFreq_LongCond, Targ_WrdFreq_ShortCond);
    P_all(i,3) = p; tvalue_all(i,3) = stats.tstat; sd_all(i,3) = stats.sd;
    
    Targ_WrdLeng_LongCond = WrdLenFreq_Targ(longid,1);
    Targ_WrdLeng_ShortCond = WrdLenFreq_Targ(shortid,1);
    [~,p,~,stats] = ttest(Targ_WrdLeng_LongCond, Targ_WrdLeng_ShortCond);
    P_all(i,4) = p; tvalue_all(i,4) = stats.tstat; sd_all(i,4) = stats.sd;
    
    PosTarg_WrdFreq_LongCond = WrdLenFreq_PosTarg(longid,3);
    PosTarg_WrdFreq_ShortCond = WrdLenFreq_PosTarg(shortid,3);
    [~,p,~,stats] = ttest(PosTarg_WrdFreq_LongCond, PosTarg_WrdFreq_ShortCond);
    P_all(i,5) = p; tvalue_all(i,5) = stats.tstat; sd_all(i,5) = stats.sd;
    
    PosTarg_WrdLeng_LongCond = WrdLenFreq_PosTarg(longid,1);
    PosTarg_WrdLeng_ShortCond = WrdLenFreq_PosTarg(shortid,1);
    [~,p,~,stats] = ttest(PosTarg_WrdLeng_LongCond, PosTarg_WrdLeng_ShortCond);
    P_all(i,6) = p; tvalue_all(i,6) = stats.tstat; sd_all(i,6) = stats.sd;
end
P_avg = mean(P_all)
stats.tstat = tvalue_all;
stats.p = P_all;
stats.sd = sd_all;
stats.p_avg = P_avg;
CondTargLeng_T_1000_loop = stats;
CondTargLeng_T_1000_loop.hdr = {'PreTarg_WrdFreq','PreTarg_WrdLeng','Targ_WrdFreq','Targ_WrdLeng','PosTarg_WrdFreq','PosTarg_WrdLeng'};
eval(['WordLength.' DataSets{dd} '.CondTargLeng_T_1000_loop = CondTargLeng_T_1000_loop;']);

%%% for single version: P_avg = [0.3605    0.7370    0.1729    0.0000    0.2981    0.4497]
%%% for combined version: P_avg = [0.2095    0.7030    0.6856    0.0000    0.2502    0.7202]
%%% for single version when loop is 30(number of subjects),P_avg=0.3687 0.7333 0.1659 0.0000 0.3087 0.4097



%% %%%======= stat for Combined datasets: get conditions based on target-length ========%%%%%
%%% target-mid-length is 6, so longer than 6 is target-long-condition,
%%% shorter than 6 is target-short-condition
load('Z:\Analyse_data\WordLength.mat');
WrdLenFreq_PreTarg = [WordLength.Targ60_1.WrdLen_PreTarg;  WordLength.JEP60.WrdLen_PreTarg];
WrdLenFreq_Targ =    [WordLength.Targ60_1.WrdLen_Targ;     WordLength.JEP60.WrdLen_Targ];
WrdLenFreq_PosTarg = [WordLength.Targ60_1.WrdLen_PosTarg;  WordLength.JEP60.WrdLen_PosTarg];
%%
TargLeng_mid = 6;
loops = 10000;
P_all = zeros(loops,6);
tvalue_all = zeros(loops,6);
sd_all = zeros(loops,6);
longid = find(WrdLenFreq_Targ(:,1)>TargLeng_mid); %% 50
shortid_all = find(WrdLenFreq_Targ(:,1)<TargLeng_mid);   %% 119
for i = 1:loops
    %     tmp_id = randperm(length(longid_all))';
    %     longid = tmp_id(1:length(shortid));
    
    tmp_id = randperm(length(shortid_all))';
    shortid = tmp_id(1:length(longid));
    
    
    PreTarg_WrdFreq_LongCond = WrdLenFreq_PreTarg(longid,3);
    PreTarg_WrdFreq_ShortCond = WrdLenFreq_PreTarg(shortid,3);
    [~,p,~,stats] = ttest(PreTarg_WrdFreq_LongCond, PreTarg_WrdFreq_ShortCond);
    P_all(i,1) = p; tvalue_all(i,1) = stats.tstat; sd_all(i,1) = stats.sd;
    
    PreTarg_WrdLeng_LongCond = WrdLenFreq_PreTarg(longid,1);
    PreTarg_WrdLeng_ShortCond = WrdLenFreq_PreTarg(shortid,1);
    [~,p,~,stats] = ttest(PreTarg_WrdLeng_LongCond, PreTarg_WrdLeng_ShortCond);
    P_all(i,2) = p; tvalue_all(i,2) = stats.tstat; sd_all(i,2) = stats.sd;
    
    Targ_WrdFreq_LongCond = WrdLenFreq_Targ(longid,3);
    Targ_WrdFreq_ShortCond = WrdLenFreq_Targ(shortid,3);
    [~,p,~,stats] = ttest(Targ_WrdFreq_LongCond, Targ_WrdFreq_ShortCond);
    P_all(i,3) = p; tvalue_all(i,3) = stats.tstat; sd_all(i,3) = stats.sd;
    
    Targ_WrdLeng_LongCond = WrdLenFreq_Targ(longid,1);
    Targ_WrdLeng_ShortCond = WrdLenFreq_Targ(shortid,1);
    [~,p,~,stats] = ttest(Targ_WrdLeng_LongCond, Targ_WrdLeng_ShortCond);
    P_all(i,4) = p; tvalue_all(i,4) = stats.tstat; sd_all(i,4) = stats.sd;
    
    PosTarg_WrdFreq_LongCond = WrdLenFreq_PosTarg(longid,3);
    PosTarg_WrdFreq_ShortCond = WrdLenFreq_PosTarg(shortid,3);
    [~,p,~,stats] = ttest(PosTarg_WrdFreq_LongCond, PosTarg_WrdFreq_ShortCond);
    P_all(i,5) = p; tvalue_all(i,5) = stats.tstat; sd_all(i,5) = stats.sd;
    
    PosTarg_WrdLeng_LongCond = WrdLenFreq_PosTarg(longid,1);
    PosTarg_WrdLeng_ShortCond = WrdLenFreq_PosTarg(shortid,1);
    [~,p,~,stats] = ttest(PosTarg_WrdLeng_LongCond, PosTarg_WrdLeng_ShortCond);
    P_all(i,6) = p; tvalue_all(i,6) = stats.tstat; sd_all(i,6) = stats.sd;
end
P_avg = mean(P_all)
stats.tstat = tvalue_all;
stats.p = P_all;
stats.sd = sd_all;
stats.p_avg = P_avg;
CondTargLeng_T_10000_loop = stats;
CondTargLeng_T_10000_loop.hdr = {'PreTarg_WrdFreq','PreTarg_WrdLeng','Targ_WrdFreq','Targ_WrdLeng','PosTarg_WrdFreq','PosTarg_WrdLeng'};
WordLength.Combine2datasets.CondTargLeng_T_1000_loop = CondTargLeng_T_1000_loop;
save('Z:\Analyse_data\WordLength','WordLength');


%% %%%================== get sentence and pre-target id based on long/short target length
load('U:\ReadingProject\stimuli\JEP_stimuli\WordsFrequencyStat_JEP.mat'); % AllWords
load('U:\ReadingProject\PTB_codes\Reading_Exp\TarWodLocFreq_JEP_1.mat'); %% same for both versions
PreTarg_id = []; %%[sentid wordid Cond_longorshort]
Targ_id = []; %%[sentid wordid Cond_longorshort]
PosTarg_id = []; %%[sentid wordid Cond_longorshort]
n = size(TarWodLocFreq,1);
Targloc = [transpose(1:n) TarWodLocFreq(:,[1 3])];
%%% condition:[long short]--[1 2]
%%% only 2 length of targets: 5 and 6 letters
for ppp = 1:n
    sid = Targloc(ppp,1);
    targloc1 = Targloc(ppp,2);
    targlength1 = AllWords.V1_length(sid,targloc1);
    cond_id = double(targlength1==5) + 1; %% [long short]--[1 2]
    PreTarg_id = [PreTarg_id; sid targloc1-1 cond_id];
    Targ_id    = [Targ_id;    sid targloc1   cond_id];
    PosTarg_id = [PosTarg_id; sid targloc1+1 cond_id];
    
    targloc2 = Targloc(ppp,3);
    targlength2 = AllWords.V1_length(sid,targloc2);
    cond_id = double(targlength2==5) + 1; %% [long short]--[1 2]
    PreTarg_id = [PreTarg_id; sid targloc2-1 cond_id];
    Targ_id    = [Targ_id;    sid targloc2   cond_id];
    PosTarg_id = [PosTarg_id; sid targloc2+1 cond_id];
end
CondTargLeng_id{1} = PreTarg_id;
CondTargLeng_id{2} = Targ_id;
CondTargLeng_id{3} = PosTarg_id;
WordLength.JEP60.CondTargLeng_id = CondTargLeng_id;
WordLength.JEP60.CondTargLeng_id_hdr = [{'PreTarg','Targ','PosTarg'};{'sentence_id','word_loc','condition_id_1_long_2_short'};];
save('Z:\Analyse_data\WordLength','WordLength');


%% %%%================== get sentence reading time
%%% using the fixation_on of the first fixation on word in the sentence and the fixation_off of the last fixation on word in the sentence as the reading time for the whole sentence
%%% note: the begining and end of the fixations should be for words in
%%% the sentence, NaN is not counted, but NaN bwteen these sentence is
%%% okay
DataSets = {'Targ60','JEP60'}; %%% now only for JEP
load('Z:\Analyse_data\SubInfo.mat');
for ddd = 1:length(DataSets)
    eval(['subs = SubInfo.subjects.' DataSets{ddd} ';']);
    eval(['PTBs = SubInfo.PTBFiles.' DataSets{ddd} ';']);
    SentDur_perwrd = zeros(length(subs),1);
    for sss = 1:length(subs)
        load(['Z:\Analyse_data\' subs{sss} '\EyeData.mat']);
        load(['Z:\RawData\PTB_data\' PTBs{sss} '.mat']);
        tmp_allfixdu = EyeData.WordFixData(:,4:3:end); %%[fix_duration]
        %%% divided by the number of fixations in each sentence, exclding
        %%% the word-skips
        wrdnum = sum(~cellfun(@isempty,Para.SentMat,'Uni',true),2);
        sentdur_tmp = nanmean(tmp_allfixdu,2);
        SentDur_perwrd(sss,1) = nanmean(sentdur_tmp);
% %         %%% divided by the number of words in each sentence
% %         wrdnum = sum(~cellfun(@isempty,Para.SentMat,'Uni',true),2);
% %         sentdur_tmp = nansum(tmp_allfixdu,2);
% %         SentDur_perwrd(sss,1) = mean(sentdur_tmp./wrdnum);
    end
    load(['Z:\Results\Behavioral\' DataSets{ddd} '_BehaData']);
    BehaData.SentDur_perwrd = SentDur_perwrd;
    save(['Z:\Results\Behavioral\' DataSets{ddd} '_BehaData'],'BehaData');
end

% % % %%%% get the averaged first-fixation-duration for words, excluding the
% % % %%%% first word.
% % % DataSets = {'Targ60','JEP60'}; %%% now only for JEP
% % % load('Z:\Analyse_data\SubInfo.mat');
% % % for ddd = 1:length(DataSets)
% % %     eval(['subs = SubInfo.subjects.' DataSets{ddd} ';']);
% % %     eval(['PTBs = SubInfo.PTBFiles.' DataSets{ddd} ';']);
% % %     SentDur_perwrd = zeros(length(subs),1);
% % %     for sss = 1:length(subs)
% % %         load(['Z:\Analyse_data\' subs{sss} '\EyeData.mat']);
% % %         load(['Z:\RawData\PTB_data\' PTBs{sss} '.mat']);
% % %         %%%% remove the first word
% % %         tmp_allfixdu = EyeData.WordFixData(:,7:3:end); %%[fix_duration]
% % %         wrdnum = sum(~cellfun(@isempty,Para.SentMat,'Uni',true),2);
% % %         sentdur_tmp = nansum(tmp_allfixdu,2);
% % %         SentDur_perwrd(sss,1) = mean(sentdur_tmp./wrdnum);
% % %     end
% % %     load(['Z:\Results\Behavioral\' DataSets{ddd} '_BehaData']);
% % %     BehaData.SentDur_perwrd_no1st = SentDur_perwrd;
% % %     save(['Z:\Results\Behavioral\' DataSets{ddd} '_BehaData'],'BehaData');
% % % end

%% get the averaged first-fixation-duration for words, excluding pre-targets
DataSets = {'Targ60','JEP60'}; %%% now only for JEP
load('Z:\Analyse_data\SubInfo.mat');
for ddd = 1:length(DataSets)
    eval(['subs = SubInfo.subjects.' DataSets{ddd} ';']);
    eval(['PTBs = SubInfo.PTBFiles.' DataSets{ddd} ';']);
    SentDur_perwrd = zeros(length(subs),1);
    for sss = 1:length(subs)
        load(['Z:\Analyse_data\' subs{sss} '\EyeData.mat']);
        load(['Z:\RawData\PTB_data\' PTBs{sss} '.mat']);
        %%%% remove pre-target word
        tmp_allfixdu = EyeData.WordFixData(:,4:3:end); %%[fix_duration]
        %%% divided by the number of fixations in each sentence, exclding
        %%% the word skipping
        if ddd == 1
           pre_id = Para.TarWodLocFreq(:,1)-1;
        else
           pre_id = Para.TarWodLocFreq(:,[1 3])-1; 
        end
        for ccc = 1:size(tmp_allfixdu,1)
            tmp_allfixdu(ccc,pre_id(ccc,:)) = nan;
        end
        sentdur_tmp = nanmean(tmp_allfixdu,2);
        SentDur_perwrd(sss,1) = nanmean(sentdur_tmp);
    end
    load(['Z:\Results\Behavioral\' DataSets{ddd} '_BehaData']);
    BehaData.SentDur_perwrd_nopre = SentDur_perwrd;
    save(['Z:\Results\Behavioral\' DataSets{ddd} '_BehaData'],'BehaData');
end




%% %==== for the cmb dataset
load('Z:\Analyse_data\SubInfo.mat');
Targ60_BD = load('Z:\Results\Behavioral\Targ60_BehaData');
JEP60_BD = load('Z:\Results\Behavioral\JEP60_BehaData');
%%%+++ get the subs that have 2 datasets
subtarg = SubInfo.subjects.Targ60;
subjep = SubInfo.subjects.JEP60;
BehaData = [];
loop_1 = {'FirstFix','Gaze','TotalGaze','PupilSize'};
loop_2 = {'Pre','Tag','Pos'};
for tmpsub = 1:length(subjep)  %%% sub id in JEP60
    stmp = find(strcmp(subjep{tmpsub}(1:end-4),subtarg)); %%% sub id in Targ60
    if stmp
        for p1 = 1:length(loop_1)
            for p2 = 1:length(loop_2)
                eval(['BehaData.' loop_1{p1} '.' loop_2{p2} '(stmp,:) = mean([Targ60_BD.BehaData.' loop_1{p1} '.' loop_2{p2} '(stmp,:); JEP60_BD.BehaData.' loop_1{p1} '.' loop_2{p2} '(tmpsub,:)],1);']);
            end
        end
         BehaData.SentDur_perwrd(stmp,1) = mean([Targ60_BD.BehaData.SentDur_perwrd(stmp,1); JEP60_BD.BehaData.SentDur_perwrd(tmpsub,1)]);
         BehaData.SentDur_perwrd_no1st(stmp,1) = mean([Targ60_BD.BehaData.SentDur_perwrd_no1st(stmp,1); JEP60_BD.BehaData.SentDur_perwrd_no1st(tmpsub,1)]);
        BehaData.SentDur_perwrd_nopre(stmp,1) = mean([Targ60_BD.BehaData.SentDur_perwrd_nopre(stmp,1); JEP60_BD.BehaData.SentDur_perwrd_nopre(tmpsub,1)]);
    end
end
%%% t-test for subjects with 2 datasets
fldnam = {'FirstFix','Gaze','TotalGaze','PupilSize'};
connam = {'Pre','Tag','Pos'};
for fff = 1:length(fldnam)
    for ccc = 1:length(connam)
        eval(['tmpdata = BehaData.' fldnam{fff} '.' connam{ccc} ';']);
        %%% remove subjects without 2 datasets
        tmpdata(tmpdata(:,1) == 0,:) = [];
        [~,p,~,stats] = ttest(tmpdata(:,1), tmpdata(:,2));
        stats.p_value = p;
        disp([fldnam{fff} '.' connam{ccc} '--------']);
        stats
        eval(['BehaData.' fldnam{fff} '.' connam{ccc} '_stat = stats;']);
    end
end
save('Z:\Results\Behavioral\CmbDatasets_BehaData','BehaData');



%%%%%%%%%%+++++++++++  acc for yes-or-no question ++++++++++%%%%%%%%%%
rootdir = 'Z:\';
load([rootdir 'Results' filesep 'Coh_hilbert_01' filesep 'minirt_trl' filesep 'Cmb.mat'])
subjects = TagCoh.subs(TagCoh.SigSubID);
load([rootdir 'Analyse_data' filesep 'SubInfo.mat']); %%% get SubInfo.MEGcode_MRIcode
PTBpath = 'Z:\RawData\PTB_data\';
acc = [];
for sss = 1:length(subjects)
   subid = find(strcmp(SubInfo.subjects.Targ60,subjects{sss}));
   ptb = SubInfo.PTBFiles.Targ60{subid};
   load([PTBpath ptb]);
   tmp = mean(Result.CORR(~isnan(Result.CORR)));
   display(['*********' subjects{sss} '---' num2str(tmp) '******']);
   acc = [acc; tmp];
   subid_jep = find(strcmp(SubInfo.subjects.JEP60,[subjects{sss} '_JEP']));
   ptb = SubInfo.PTBFiles.JEP60{subid_jep};
   load([PTBpath ptb]);
   tmp = mean(Result.CORR(~isnan(Result.CORR)));
   display(['*********JEP ' subjects{sss} '---' num2str(tmp) '******']);
   acc = [acc; tmp];
end
meanacc = mean(acc);
stdacc = std(acc);
load('Z:\Results\Behavioral\CmbDatasets\CmbDatasets_BehaData')
BehaData.YesNoQ = [meanacc stdacc];
save('Z:\Results\Behavioral\CmbDatasets\CmbDatasets_BehaData','BehaData')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% bi/trigram freq, orth neighb %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get all the targets
ppath = 'Z:\Analyse_codes\PTB_codes\Reading_Exp\';
%======Targ60 (every sub read all targets, so stimuli_version doesn't matter)
load([ppath 'SentMat_1.mat']) % SentMat
load([ppath 'TarWodLocFreq_1.mat']) % TarWodLocFreq
n = size(TarWodLocFreq,1);
words = cell(n,1); % [low high]condition
for s = 1:n
    tgid = TarWodLocFreq(s,1); % target_id
    words(s,:) = SentMat(s, tgid);
end
words_low_high = [words(TarWodLocFreq(:,2)==1) words(TarWodLocFreq(:,2)==2)];
TargetVariables.Targ60.words_low_high = words_low_high;
%%% copy words into N_watch, then copy values back to [low high]. 
%%% Note, N-watch will lost the last word when 'copy-paste'
words_low_high(end+1,:) = [{''},{''}];
vbs = {'CELEX','BFTK','TFTK','NeigbSiz'};
for vv = 1:length(vbs)
    [~,p,~,stats] = ttest(low(:,vv),high(:,vv));
    stats.p_value = p;
    disp([vbs{vv} ': p is ' num2str(p)]);
    eval(['TargetVariables.Targ60.' vbs{vv} '_stats = stats;']);
end
TargetVariables.Targ60.low_CELEX_BF_TF_NS = low;
TargetVariables.Targ60.high_CELEX_BF_TF_NS = high;


%===== JEP60 (two sentence_versions)
% version_1
load([ppath 'SentMat_JEP_1.mat']) % SentMat
load([ppath 'TarWodLocFreq_JEP_1.mat']) % TarWodLocFreq
n = size(TarWodLocFreq,1);
%%% note: two words in one sentence
TarWodLocFreq = [TarWodLocFreq; TarWodLocFreq];
TarWodLocFreq(n+1:end,1) = TarWodLocFreq(1:n,3);
TarWodLocFreq(:,3) = repmat([1:n]',2,1);
n = size(TarWodLocFreq,1);
words = cell(n,1); % [low high]condition
for s = 1:n
    tgid = TarWodLocFreq(s,1); % target_id
    sentid = TarWodLocFreq(s,3); % sentence_id
    words(s,:) = SentMat(sentid, tgid);
end
words_low_high = [words(TarWodLocFreq(:,2)==1) words(TarWodLocFreq(:,2)==2)];
% version_2
load([ppath 'SentMat_JEP_2.mat']) % SentMat
load([ppath 'TarWodLocFreq_JEP_2.mat']) % TarWodLocFreq
n = size(TarWodLocFreq,1);
TarWodLocFreq = [TarWodLocFreq; TarWodLocFreq];
TarWodLocFreq(n+1:end,1) = TarWodLocFreq(1:n,3);
TarWodLocFreq(:,3) = repmat([1:n]',2,1);
n = size(TarWodLocFreq,1);
words = cell(n,1); % [low high]condition
for s = 1:n
    tgid = TarWodLocFreq(s,1); % target_id
    sentid = TarWodLocFreq(s,3); % sentence_id
    words(s,:) = SentMat(sentid, tgid);
end
words_low_high_2 = [words(TarWodLocFreq(:,2)==1) words(TarWodLocFreq(:,2)==2)];
% combine version_1 and 2
words_low_high = [words_low_high; words_low_high_2];
TargetVariables.JEP60.words_low_high = words_low_high;
%%% copy words into N_watch, then copy values back to [low high]. 
%%% Note, N-watch will lost the last word when 'copy-paste'
words_low_high(end+1,:) = [{''},{''}];
vbs = {'CELEX','BFTK','TFTK','NeigbSiz'};
for vv = 1:length(vbs)
    [~,p,~,stats] = ttest(low(:,vv),high(:,vv));
    stats.p_value = p;
    disp([vbs{vv} ': p is ' num2str(p)]);
    eval(['TargetVariables.JEP60.' vbs{vv} '_stats = stats;']);
end
TargetVariables.JEP60.low_CELEX_BF_TF_NS = low;
TargetVariables.JEP60.high_CELEX_BF_TF_NS = high;


%%%======= combine Targ60 and JEP60
words_low_high = [TargetVariables.Targ60.words_low_high; TargetVariables.JEP60.words_low_high];
TargetVariables.combined.words_low_high = words_low_high;
words_low_high(end+1,:) = [{''},{''}];
%%% copy words into N-watch to get measurements, which are copied back into
low = [];
high = [];
vbs = {'CELEX','BFtoken','TFtoken','NeigbSiz','BFtype','TFtype'};
for vv = 1:length(vbs)
    [~,p,~,stats] = ttest(low(:,vv),high(:,vv));
    stats.p_value = p;
    disp([vbs{vv} ': p is ' num2str(p)]);
    eval(['TargetVariables.combined.' vbs{vv} '_stats = stats;']);
end
TargetVariables.combined.low_CELEX_BFtoken_TFtoken_NS_BFtype_TFtype = low;
TargetVariables.combined.high_CELEX_BFtoken_TFtoken_NS_BFtype_TFtype = high;

% correlation between bi/trigram and word freq
all_var = [low; high];
% word_freq & bigram_freq --- [0.2885 9.1166e-11]
[coef, pval] = corr(all_var(:,1), all_var(:,2),'type','Pearson')
TargetVariables.combined.CorrP_CELEX_BF = [coef, pval];
% word_freq & trigram_freq --- [0.4413 1.4061e-24]
[coef, pval] = corr(all_var(:,1), all_var(:,3),'type','Pearson')
TargetVariables.combined.CorrP_CELEX_TF = [coef, pval];
% bigram_freq & trigram_freq --- [0.6755 4.8845e-66]
[coef, pval] = corr(all_var(:,2), all_var(:,3),'type','Pearson')
TargetVariables.combined.CorrP_BF_TF = [coef, pval];
% word_freq & NeigbSiz --- [0.0576 0.2048]
[coef, pval] = corr(all_var(:,1), all_var(:,4),'type','Pearson')
TargetVariables.combined.CorrP_CELEX_NS = [coef, pval];


need to get the bi/trigram freq for the pre-targets in Targ60 as well!

save('Z:\Results\Behavioral\TargetVariables','TargetVariables');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% bi/trigram freq, orth neighb_seperate by version %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get all the targets
ppath = 'Z:\Analyse_codes\PTB_codes\Reading_Exp\';
% loop over version
%%% Targ60
for v = [1 2]
    load([ppath 'SentMat_' num2str(v) '.mat']) % SentMat
    load([ppath 'TarWodLocFreq_' num2str(v) '.mat']) % TarWodLocFreq
    n = size(TarWodLocFreq,1);
    words = cell(n,1); % [low high]condition
    for s = 1:n
        tgid = TarWodLocFreq(s,1); % target_id
        words(s,:) = SentMat(s, tgid);
    end
    eval(['TargetVariables.Targ60_v' num2str(v) '.words = words;']);
end
%%% JEP (note: two target words per sentence)
for v = [1 2]
    load([ppath 'SentMat_JEP_' num2str(v) '.mat']) % SentMat
    load([ppath 'TarWodLocFreq_JEP_' num2str(v) '.mat']) % TarWodLocFreq
    n = size(TarWodLocFreq,1);
    TarWodLocFreq = [TarWodLocFreq; TarWodLocFreq];
    TarWodLocFreq(n+1:end,1) = TarWodLocFreq(1:n,3);
    TarWodLocFreq(:,3) = repmat([1:n]',2,1);
    n = size(TarWodLocFreq,1);
    words = cell(n,1); % [low high]condition
    for s = 1:n
        tgid = TarWodLocFreq(s,1); % target_id
        sentid = TarWodLocFreq(s,3);
        words(s,:) = SentMat(sentid, tgid);
    end
    eval(['TargetVariables.JEP60_v' num2str(v) '.words = words;']);
    eval(['TargetVariables.JEP60_v' num2str(v) '.SentID_TargID = TarWodLocFreq(:,[3 1]);']);
end
%%% copy words into N_watch, then copy values back 
TargetVariables.Targ60_v1.CELEX_BFtoken_TFtoken_Length_NS = ;
TargetVariables.Targ60_v2.CELEX_BFtoken_TFtoken_Length_NS = ;
TargetVariables.JEP60_v1.CELEX_BFtoken_TFtoken_Length_NS = ;
TargetVariables.JEP60_v2.CELEX_BFtoken_TFtoken_Length_NS = ;

%%% give condition_index based on ortho var
vbs = {'BFtoken','TFtoken','NeigbSiz'};
lowtmp = TargetVariables.combined.low_CELEX_BFtoken_TFtoken_NS_BFtype_TFtype_Length;
hightmp = TargetVariables.combined.high_CELEX_BFtoken_TFtoken_NS_BFtype_TFtype_Length;
alltmp = [lowtmp; hightmp];
tmp = alltmp(:,[2 3 4]);
wrdlen = alltmp(:,7);
CondId = nan(size(tmp));
hn = size(tmp,1)/2; % half n
for vb = 1:length(vbs)
    [~, I] = sort(tmp(:,vb), 'ascend'); %[1-low;2-high]
    CondId(I(1:hn),vb) = 1;
    CondId(I(hn+1:end),vb) = 2;
    %%% get word length ttest for each orth_type
    len_1 = wrdlen(I(1:hn));
    len_2 = wrdlen(I(hn+1:end));
    [~,p,~,stats] = ttest(len_1, len_2); p
    stats.p = p;
    eval(['TargetVariables.combined.Ttest_Length_dividedby_' vbs{vb} '= stats;']);
end
save('Z:\Results\Behavioral\TargetVariables','TargetVariables');


%%%%% select high/low trigram token frequency conditions, then run stat for
%%%%% other variables
a = [TargetVariables.combined.low_CELEX_BFtoken_TFtoken_NS_BFtype_TFtype_Length;...
     TargetVariables.combined.high_CELEX_BFtoken_TFtoken_NS_BFtype_TFtype_Length];

%%% bigram token frequency
[bitokem, I] = sort(a(:,2));
bi_low = a(I(1:243),:);
bi_high = a(I(244:end),:);
%%%
mean(bi_low)  %[29.27	520.51	97.39	1.10	34.75	6.30	5.91]
mean(bi_high) %[71.31	1791.45	376.92	3.16	80.33	14.68	5.65]

[~,p,~,stats] = ttest(bi_low(:,3), bi_high(:,3)); disp(['p_TFtoken:' num2str(p)]) %p=0
[~,p,~,stats] = ttest(bi_low(:,1), bi_high(:,1)); disp(['p_CELEX:' num2str(p)])   %p=2.2915e-05
[~,p,~,stats] = ttest(bi_low(:,4), bi_high(:,4)); disp(['p_NS:' num2str(p)])      %p=0
[~,p,~,stats] = ttest(bi_low(:,7), bi_high(:,7)); disp(['p_Length:' num2str(p)])  %p=0.0007

 
%%% trigram token frequency
[tritokem, I] = sort(a(:,3));
tri_low = a(I(1:243),:);
tri_high = a(I(244:end),:);
%%%
mean(tri_low)  %[14.76	713.06  59.05	1.28	41.74	5.71	5.77]
mean(tri_high) %[85.83 1598.89	415.27	2.97	73.34	15.27	5.78]

[~,p,~,stats] = ttest(tri_low(:,3), tri_high(:,3)); disp(['p_TFtoken:' num2str(p)]) %p=0
[~,p,~,stats] = ttest(tri_low(:,1), tri_high(:,1)); disp(['p_CELEX:' num2str(p)])   %p=2.9199e-14
[~,p,~,stats] = ttest(tri_low(:,4), tri_high(:,4)); disp(['p_NS:' num2str(p)])      %p=1.6764e-14
[~,p,~,stats] = ttest(tri_low(:,7), tri_high(:,7)); disp(['p_Length:' num2str(p)])  %p=0.90925












%%% analyse of behavioral the eye movement data
%%% edit: YPan, 20220126

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%========== stat on behavioral data ========== %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ddd = 1; % index for task

%%% paths
PPath.DataPath = 'Z:\Semantic\Analyse_data\';
load([PPath.DataPath 'ExpInfo.mat']) % ExpInfo
DS = ExpInfo.DSName{ddd};
eval(['subjects = ExpInfo.subjects.' DS ';']);
PPath.FigPath = ['Z:\Semantic\Results\' DS '\Behav\']; %result path
if (~exist(PPath.FigPath,'dir'))
    mkdir(PPath.FigPath);
end
addpath(genpath('Z:\Semantic\Analyse_codes\plotfig\'))

%%% remove empty row in subjects
subjects = subjects(~strcmp(subjects,''));
CondId = ExpInfo.CondID{ddd};

%% % gaze duration: all fixation on target before eye move out of taget,
%%% re-read is not included
BehaData = [];
BehaData.CondName = ExpInfo.CondName{ddd};
for sss = 1:length(subjects)
    load([PPath.DataPath DS '_' subjects{sss} '\Event.mat']); %Event
    if sss == 1
        % get the column id for some variable
        fixdur = find(strcmp(Event.event_raw_header,'fixation_duration'));
        loc2targ = find(strcmp(Event.event_raw_header,'loc2targ'));
        firstPass = find(strcmp(Event.event_raw_header,'FirstPassFix'));
        cond = find(strcmp(Event.event_raw_header,'SentenceCondition'));
        pupil = find(strcmp(Event.event_raw_header,'PupilSize'));
        preorder = find(strcmp(Event.event_raw_header,'PreviousOrder'));
    end
    
    tmp = Event.event_raw;
    %%{'sentence_id'},{'word_loc'},{'loc2targ'},{'word_freq'},{'word_length'},{'saccade2this_du…'},{'fixation_on_MEG'},{'fixation_duration'},{'NextOrder'},{'FirstPassFix'},{'PreviousOrder'},{'SentenceCond'},{'PupilSize'}
    %%% remove outliers
    tmp = tmp(tmp(:,fixdur)>=80 & tmp(:,fixdur)<=1000,:);
    
    %%% get duration & pupil size: [low high]
    % index for first fixation
    ff = tmp(:,firstPass)==1;
    pre1 = tmp(:,loc2targ)==-1 & tmp(:,cond)==CondId(1);
    pre1_ff = tmp(pre1 & ff,[fixdur pupil]);
    pre2 = tmp(:,loc2targ)==-1 & tmp(:,cond)==CondId(2);
    pre2_ff = tmp(pre2 & ff,[fixdur pupil]);
    tag1 = tmp(:,loc2targ)==0 & tmp(:,cond)==CondId(1);
    tag1_ff = tmp(tag1 & ff,[fixdur pupil]);
    tag2 = tmp(:,loc2targ)==0 & tmp(:,cond)==CondId(2);
    tag2_ff = tmp(tag2 & ff,[fixdur pupil]);
    pos1 = tmp(:,loc2targ)==1 & tmp(:,cond)==CondId(1);
    pos1_ff = tmp(pos1 & ff,[fixdur pupil]);
    pos2 = tmp(:,loc2targ)==1 & tmp(:,cond)==CondId(2);
    pos2_ff = tmp(pos2 & ff,[fixdur pupil]);
    %%% duration of the first fixation 
    FirstFix.Pre(sss,:) = [mean(pre1_ff(:,1))   mean(pre2_ff(:,1))];
    FirstFix.Tag(sss,:) = [mean(tag1_ff(:,1))   mean(tag2_ff(:,1))];
    FirstFix.Pos(sss,:) = [mean(pos1_ff(:,1))   mean(pos2_ff(:,1))];
    %%% number of the first fixation 
    FirstFix.Pre_Num(sss,:) = [size(pre1_ff,1)   size(pre2_ff,1)];
    FirstFix.Tag_Num(sss,:) = [size(tag1_ff,1)   size(tag2_ff,1)];
    FirstFix.Pos_Num(sss,:) = [size(pos1_ff,1)   size(pos2_ff,1)];
    %%% pupil size of the first fixation
    PupilSize.Pre(sss,:) = [mean(pre1_ff(:,2))   mean(pre2_ff(:,2))];
    PupilSize.Tag(sss,:) = [mean(tag1_ff(:,2))   mean(tag2_ff(:,2))];
    PupilSize.Pos(sss,:) = [mean(pos1_ff(:,2))   mean(pos2_ff(:,2))];
    
    %%% gaze: sum of the fixations before eyes leave a given word, the saccades duration is not included here   
    gaz = ff | tmp(:,preorder)==0;
    Gaze.Pre(sss,:) = [sum(tmp(pre1&gaz,fixdur))./size(pre1_ff,1),...
                       sum(tmp(pre2&gaz,fixdur))./size(pre2_ff,1)];
    Gaze.Tag(sss,:) = [sum(tmp(tag1&gaz,fixdur))./size(tag1_ff,1),...
                       sum(tmp(tag2&gaz,fixdur))./size(tag2_ff,1)];
    Gaze.Pos(sss,:) = [sum(tmp(pos1&gaz,fixdur))./size(pos1_ff,1),...
                       sum(tmp(pos2&gaz,fixdur))./size(pos2_ff,1)];
    %%% number of gaze
    Gaze.Pre_Num(sss,:) = [sum(pre1&gaz == 1) sum(pre2&gaz == 1)];
    Gaze.Tag_Num(sss,:) = [sum(tag1&gaz == 1) sum(tag2&gaz == 1)];
    Gaze.Pos_Num(sss,:) = [sum(pos1&gaz == 1) sum(pos2&gaz == 1)];
           
    %%% total gaze: sum of all the fixations to a given word, the saccades duration is not included here  
    TotalGaze.Pre(sss,:) = [sum(tmp(pre1,fixdur))./size(pre1_ff,1),...
                       sum(tmp(pre2,fixdur))./size(pre2_ff,1)];
    TotalGaze.Tag(sss,:) = [sum(tmp(tag1,fixdur))./size(tag1_ff,1),...
                       sum(tmp(tag2,fixdur))./size(tag2_ff,1)];
    TotalGaze.Pos(sss,:) = [sum(tmp(pos1,fixdur))./size(pos1_ff,1),...
                       sum(tmp(pos2,fixdur))./size(pos2_ff,1)];
    %%% number of total gaze
    TotalGaze.Pre_Num(sss,:) = [sum(pre1 == 1) sum(pre2 == 1)];
    TotalGaze.Tag_Num(sss,:) = [sum(tag1 == 1) sum(tag2 == 1)];
    TotalGaze.Pos_Num(sss,:) = [sum(pos1 == 1) sum(pos2 == 1)];
end

%% simple paired t-tests
% durations of first fixations
[~,p,~,stat] = ttest(FirstFix.Pre(:,1),FirstFix.Pre(:,2));
stat.p = p
FirstFix.Pre_stat = stat;
[~,p,~,stat] = ttest(FirstFix.Tag(:,1),FirstFix.Tag(:,2));
stat.p = p
FirstFix.Tag_stat = stat;
[~,p,~,stat] = ttest(FirstFix.Pos(:,1),FirstFix.Pos(:,2));
stat.p = p
FirstFix.Pos_stat = stat;
% durations of gaze fixations
[~,p,~,stat] = ttest(Gaze.Pre(:,1),Gaze.Pre(:,2));
stat.p = p
Gaze.Pre_stat = stat;
[~,p,~,stat] = ttest(Gaze.Tag(:,1),Gaze.Tag(:,2));
stat.p = p
Gaze.Tag_stat = stat;
[~,p,~,stat] = ttest(Gaze.Pos(:,1),Gaze.Pos(:,2));
stat.p = p
Gaze.Pos_stat = stat;
% durations of total gaze 
[~,p,~,stat] = ttest(TotalGaze.Pre(:,1),TotalGaze.Pre(:,2));
stat.p = p
TotalGaze.Pre_stat = stat;
[~,p,~,stat] = ttest(TotalGaze.Tag(:,1),TotalGaze.Tag(:,2));
stat.p = p
TotalGaze.Tag_stat = stat;
[~,p,~,stat] = ttest(TotalGaze.Pos(:,1),TotalGaze.Pos(:,2));
stat.p = p
TotalGaze.Pos_stat = stat;
% pupil size of first fixation
[~,p,~,stat] = ttest(PupilSize.Pre(:,1),PupilSize.Pre(:,2));
stat.p = p
PupilSize.Pre_stat = stat;
[~,p,~,stat] = ttest(PupilSize.Tag(:,1),PupilSize.Tag(:,2));
stat.p = p
PupilSize.Tag_stat = stat;
[~,p,~,stat] = ttest(PupilSize.Pos(:,1),PupilSize.Pos(:,2));
stat.p = p
PupilSize.Pos_stat = stat;

% percentage of non-first-fixation gaze. Note: number of first fixation
% indicates the number of fixated words, which used as denominator here
n_nonff = (Gaze.Pre_Num - FirstFix.Pre_Num)./ FirstFix.Pre_Num;
[~,p,~,stat] = ttest(n_nonff(:,1),n_nonff(:,2));
stat.p = p
Gaze.Pre_nonfirst_percent = n_nonff;
Gaze.Pre_nonfirst_percent_stat = stat;
n_nonff = (Gaze.Tag_Num - FirstFix.Tag_Num)./ FirstFix.Tag_Num;
[~,p,~,stat] = ttest(n_nonff(:,1),n_nonff(:,2));
stat.p = p
Gaze.Tag_nonfirst_percent = n_nonff;
Gaze.Tag_nonfirst_percent_stat = stat;
n_nonff = (Gaze.Pos_Num - FirstFix.Pos_Num)./ FirstFix.Pos_Num;
[~,p,~,stat] = ttest(n_nonff(:,1),n_nonff(:,2));
stat.p = p
Gaze.Pos_nonfirst_percent = n_nonff;
Gaze.Pos_nonfirst_percent_stat = stat;

% percentage of regression in total gaze; backward eye-movements to earlier points in the text
n_reg = (TotalGaze.Pre_Num - Gaze.Pre_Num)./ FirstFix.Pre_Num;
[~,p,~,stat] = ttest(n_reg(:,1),n_reg(:,2));
stat.p = p
Regression_percent.Pre = n_reg;
Regression_percent.Pre_stat = stat;
n_reg = (TotalGaze.Tag_Num - Gaze.Tag_Num)./ FirstFix.Tag_Num;
[~,p,~,stat] = ttest(n_reg(:,1),n_reg(:,2));
stat.p = p
Regression_percent.Tag = n_reg;
Regression_percent.Tag_stat = stat;
n_reg = (TotalGaze.Pos_Num - Gaze.Pos_Num)./ FirstFix.Pos_Num;
[~,p,~,stat] = ttest(n_reg(:,1),n_reg(:,2));
stat.p = p
Regression_percent.Pos = n_reg;
Regression_percent.Pos_stat = stat;

%% save out
BehaData.FirstFix = FirstFix;
BehaData.Gaze = Gaze;
BehaData.TotalGaze = TotalGaze;
BehaData.PupilSize = PupilSize;
BehaData.Regression_percent = Regression_percent;
save([PPath.FigPath 'BehaData.mat'],'BehaData');

%% %%%%%%%%%%%%%=============== plotting =================%%%%%%%%%%%%%%%%
FigNam = {'FirstFix';'Gaze';'TotalGaze';'PupilSize';'Regression_percent'};
for ff = 1:length(FigNam)
    eval(['tmp = BehaData.' FigNam{ff} ';']);    
    h = figure('Name',FigNam{ff},'color',[1 1 1]);
    y = [mean(tmp.Pre); mean(tmp.Tag); mean(tmp.Pos)];         % random y values (3 groups of 2 parameters)
    errY = [std(tmp.Pre); std(tmp.Tag); std(tmp.Pos)]./sqrt(size(tmp.Pre,1));
    a = barwitherr(cat(3,-errY,errY), y);    % Plot with errorbars
    if ff == 1
        ylim([160 260])
        yyy = 245;
    elseif ff == 2
        ylim([220 320])
        yyy = 295;
    elseif ff == 3
        ylim([200 450])
        yyy = 400;
    elseif ff == 4
        ylim([4600 5600])
        yyy = 5460;
    else  
        ylim([0 0.7])
        yyy = 0.57;
    end
    set(gca,'box','off','LineWidth',2)
    set(gca,'XTickLabel',{'Pre-target','Target','Post-target'},'FontSize',12,'FontWeight','bold')
    text(1,yyy,['p = ' num2str(round(1000*tmp.Pre_stat.p)/1000)],'FontSize',10,'FontWeight','bold');
    text(2,yyy,['p = ' num2str(round(1000*tmp.Tag_stat.p)/1000)],'FontSize',10,'FontWeight','bold');
    text(3,yyy,['p = ' num2str(round(1000*tmp.Pos_stat.p)/1000)],'FontSize',10,'FontWeight','bold');
    switch ff
        case 1
            subtitle = 'First Fixation duration(ms)';
        case 2
            subtitle = 'Gaze duration(ms)';
        case 3
            subtitle = 'Total viewing time(ms)';
        case 4
            subtitle = 'Pupil Dilation';
        case 5
            subtitle = 'Percent of regression';
    end
    title(subtitle,'FontSize',18,'FontWeight','bold');
    legendflex(a,[{'Incong'},{'Cong'}],'anchor', {'ne','ne'}, 'buffer', [5 -5],'fontsize',10,'FontWeight','bold','xscale',0.8,'box', 'off');
    saveas(h,[PPath.FigPath FigNam{ff}]);
end



%% %%%%%%%% acc for yes-or-no question %%%%%%%%%%
PTBpath = 'Z:\Semantic\RawData\PTB_data\';
eval(['ptb_file = ExpInfo.PTBFiles.' DS ';']);
%%% remove empty row
ptb_file = subjects(~strcmp(ptb_file,''));

acc = [];
for sss = 1:length(ptb_file)
   load([PTBpath ptb_file{sss}]);
   %%% select the trls based on task
   if isfield(Para,'TargCondPair')
       thistask = arrayfun(@(x) ismember(x,CondId),Para.TargCondPair(:,2),'Uni',true);
   else % for sub1-4
       thistask = arrayfun(@(x) ismember(x,CondId),Para.TargLoc_Cond(:,3),'Uni',true);
   end
   tmp = nanmean(Result.CORR(thistask));
   acc = [acc; tmp];
end
meanacc = mean(acc);
stdacc = std(acc);
BehaData.YesNoQ_AvgStd = [meanacc stdacc];
save([PPath.FigPath 'BehaData.mat'],'BehaData');











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
%%
BehaData.CondName = ExpInfo.CondName{ddd};
for sss = 1:length(subjects)
    load([PPath.DataPath DS '_' subjects{sss} '\Event.mat']); %Event
    if sss == 1
        % get the column id for some variable
        sentid = find(strcmp(Event.event_raw_header,'sentence_id'));
        wordloc = find(strcmp(Event.event_raw_header,'word_loc'));
        loc2targ = find(strcmp(Event.event_raw_header,'loc2targ'));
        sacdur = find(strcmp(Event.event_raw_header,'saccade2this_duration')); 
        megtrig = find(strcmp(Event.event_raw_header,'fixation_on_MEG')); 
        fixdur = find(strcmp(Event.event_raw_header,'fixation_duration'));
        firstPass = find(strcmp(Event.event_raw_header,'FirstPassFix'));%note!first pass more means the first fixation 
        cond = find(strcmp(Event.event_raw_header,'SentenceCondition'));
        pupil = find(strcmp(Event.event_raw_header,'PupilSize'));
        preorder = find(strcmp(Event.event_raw_header,'PreviousOrder'));
        %%% pre-config the gaze event
        GazeEvent.hdr = Event.event_raw_header([sentid wordloc loc2targ megtrig fixdur cond]);
    end
    
    
    %%%%%%%%%%%%%========== get gaze_event and save it out ========== %%%%%%%%%%%%%%
    %%% gaze: sum of the fixations before eyes leave a given word, the saccades duration ARE included here   
    %step1: find all the non-first fixations that in the same word and add
    %the saccade duration before it to the fixation
    tmp = Event.event_raw;
    tmp = tmp(tmp(:,fixdur)>=80,:); % get rid of the gaze that shorter than 80ms
    samewrdfix = tmp(:,preorder)==0; %index for the non-first-fixation before eyes leave this word
    tmp(samewrdfix,fixdur) = nansum([tmp(samewrdfix,fixdur) tmp(samewrdfix,sacdur)],2);
    id_firstfirx = find(tmp(:,firstPass)==1);
    n = length(id_firstfirx);
    GazeEvent.event = zeros(n,length(GazeEvent.hdr));
    for i = 1:n
        cur_i = id_firstfirx(i); %row id in the raw tmp array
        tmpgz = tmp(cur_i,fixdur);
        loopi = cur_i + 1;
        while loopi < n && tmp(loopi,preorder) == 0
            tmpgz = [tmpgz tmp(loopi,fixdur)];
            loopi = loopi + 1;
        end
        GazeEvent.event(i,:) = tmp(id_firstfirx(i),[sentid wordloc loc2targ megtrig fixdur cond]);
        GazeEvent.event(i,end-1) = nansum(tmpgz);
    end
    %%% save out
    save([PPath.DataPath DS '_' subjects{sss} '\GazeEvent.mat'],'GazeEvent');
    
    %%% get gaze data
    gz_pre1 = GazeEvent.event(:,3)==-1 & GazeEvent.event(:,6)==CondId(1);
    gz_pre2 = GazeEvent.event(:,3)==-1 & GazeEvent.event(:,6)==CondId(2);
    Gaze.Pre(sss,:) = [mean(GazeEvent.event(gz_pre1,5)) mean(GazeEvent.event(gz_pre2,5))];
    gz_tag1 = GazeEvent.event(:,3)==0 & GazeEvent.event(:,6)==CondId(1);
    gz_tag2 = GazeEvent.event(:,3)==0 & GazeEvent.event(:,6)==CondId(2);
    Gaze.Tag(sss,:) = [mean(GazeEvent.event(gz_tag1,5)) mean(GazeEvent.event(gz_tag2,5))];
    gz_pos1 = GazeEvent.event(:,3)==1 & GazeEvent.event(:,6)==CondId(1);
    gz_pos2 = GazeEvent.event(:,3)==1 & GazeEvent.event(:,6)==CondId(2);
    Gaze.Pos(sss,:) = [mean(GazeEvent.event(gz_pos1,5)) mean(GazeEvent.event(gz_pos2,5))];
    
    %%%%%%%%%%%%%========== get total gaze duration and regression probability========== %%%%%%%%%%%%%%
    %%% total gze: sum of all the fixations to a given word, the saccades duration IS included here
    togz_reg = zeros(n,4); %[loc2targ cond fixdur num_of_regression]
    for i = 1:n
        cur_i = id_firstfirx(i); %row id in the raw tmp array
        si = tmp(cur_i,sentid);
        wi = tmp(cur_i,loc2targ);
        cur_wrd = tmp(:,sentid)==si & tmp(:,loc2targ)==wi;
        togz_reg(i,1:3) = [wi tmp(cur_i,cond) sum(tmp(cur_wrd,fixdur))]; 
        %%% get the probablity of regressed fixations
        % We calculated the probability of making a regression into a word 
        % as the proportion of trials on which there was at least one regression 
        % from a later part of the sentence back to that word. 
        togz_reg(i,4) = any(tmp(cur_wrd,preorder) > 0);
    end
    % get total gaze durations in each condition
    tg_pre1 = togz_reg(:,1)==-1 & togz_reg(:,2)==CondId(1);
    tg_pre2 = togz_reg(:,1)==-1 & togz_reg(:,2)==CondId(2);
    Totalgze.Pre(sss,:) = [mean(togz_reg(tg_pre1,3)) mean(togz_reg(tg_pre2,3))];
    tg_tag1 = togz_reg(:,1)==0 & togz_reg(:,2)==CondId(1);
    tg_tag2 = togz_reg(:,1)==0 & togz_reg(:,2)==CondId(2);
    Totalgze.Tag(sss,:) = [mean(togz_reg(tg_tag1,3)) mean(togz_reg(tg_tag2,3))];
    tg_pos1 = togz_reg(:,1)==1 & togz_reg(:,2)==CondId(1);
    tg_pos2 = togz_reg(:,1)==1 & togz_reg(:,2)==CondId(2);
    Totalgze.Pos(sss,:) = [mean(togz_reg(tg_pos1,3)) mean(togz_reg(tg_pos2,3))];
    % get regression probability in each condition
    RegsProb.Pre(sss,:) = [mean(togz_reg(tg_pre1,4)) mean(togz_reg(tg_pre2,4))];
    RegsProb.Tag(sss,:) = [mean(togz_reg(tg_tag1,4)) mean(togz_reg(tg_tag2,4))];
    RegsProb.Pos(sss,:) = [mean(togz_reg(tg_pos1,4)) mean(togz_reg(tg_pos2,4))];

   
    %%%%%%%%%%%%%========== get first fixation ========== %%%%%%%%%%%%%%
    %%% remove outliers, only for the first fixation
    tmp = Event.event_raw;
    tmp = tmp(tmp(:,fixdur)>=80 & tmp(:,fixdur)<=1000,:);
    %%% get duration & pupil size: [cond1 cond2]
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
    %%% pupil size of the first fixation
    PupilSize.Pre(sss,:) = [mean(pre1_ff(:,2))   mean(pre2_ff(:,2))];
    PupilSize.Tag(sss,:) = [mean(tag1_ff(:,2))   mean(tag2_ff(:,2))];
    PupilSize.Pos(sss,:) = [mean(pos1_ff(:,2))   mean(pos2_ff(:,2))];
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
% durations of gze fixations
[~,p,~,stat] = ttest(Gaze.Pre(:,1),Gaze.Pre(:,2));
stat.p = p
Gaze.Pre_stat = stat;
[~,p,~,stat] = ttest(Gaze.Tag(:,1),Gaze.Tag(:,2));
stat.p = p
Gaze.Tag_stat = stat;
[~,p,~,stat] = ttest(Gaze.Pos(:,1),Gaze.Pos(:,2));
stat.p = p
Gaze.Pos_stat = stat;
% durations of total gze 
[~,p,~,stat] = ttest(Totalgze.Pre(:,1),Totalgze.Pre(:,2));
stat.p = p
Totalgze.Pre_stat = stat;
[~,p,~,stat] = ttest(Totalgze.Tag(:,1),Totalgze.Tag(:,2));
stat.p = p
Totalgze.Tag_stat = stat;
[~,p,~,stat] = ttest(Totalgze.Pos(:,1),Totalgze.Pos(:,2));
stat.p = p
Totalgze.Pos_stat = stat;
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

% probability of regression  
[~,p,~,stat] = ttest(RegsProb.Pre(:,1),RegsProb.Pre(:,2));
stat.p = p
RegsProb.Pre_stat = stat;
[~,p,~,stat] = ttest(RegsProb.Tag(:,1),RegsProb.Tag(:,2));
stat.p = p
RegsProb.Tag_stat = stat;
[~,p,~,stat] = ttest(RegsProb.Pos(:,1),RegsProb.Pos(:,2));
stat.p = p
RegsProb.Pos_stat = stat;


%% save out
BehaData.FirstFix = FirstFix;
BehaData.Gaze = Gaze;
BehaData.Totalgze = Totalgze;
BehaData.PupilSize = PupilSize;
BehaData.RegsProb = RegsProb;
save([PPath.FigPath 'BehaData.mat'],'BehaData');

%% %%%%%%%%%%%%%=============== plotting =================%%%%%%%%%%%%%%%%
FigNam = {'FirstFix';'Gaze';'Totalgze';'PupilSize';'RegsProb'};
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
            subtitle = 'gze duration(ms)';
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











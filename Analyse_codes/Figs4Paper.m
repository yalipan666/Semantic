%% plot figures for the manuscript of semantic parafoveal project
%  edit by Yali, 20220222

figpath = 'U:\writing\Parafoveal_Semantic\figs\';
addpath(genpath('Z:\Semantic\Analyse_codes\plotfig\'))
%%% for colormap
cmap = colormap(cbrewer('div','RdBu',32));
cmap = colormap(flipud(cmap));

%% eye movement behavioral data
cd Z:\Semantic\Results\sv\Behav\
load('BehaData.mat')
FigNam = {'FirstFix';'TotalGaze';'Regression_percent'};
for ff = 1:length(FigNam)
    eval(['tmp = BehaData.' FigNam{ff} ';']);
    if ff == 3
        tmp.Pre = tmp.Pre*100;
        tmp.Tag = tmp.Tag*100;
    end
    h = figure('Name',FigNam{ff},'color',[1 1 1],'position',[300 300 120 140]);
    y = [mean(tmp.Pre); mean(tmp.Tag);];         % random y values (2 groups of 2 parameters)
    errY = [std(tmp.Pre); std(tmp.Tag);]./sqrt(size(tmp.Pre,1));
    a = barwitherr(cat(3,-errY,errY), y);    % Plot with errorbars
    if ff == 1
        subtitle = 'First fixation duration(ms)';
        ylim([180 280])
        yyy = 245;
    elseif ff == 2
        subtitle = 'Total viewing time(ms)';
        ylim([200 450])
        yyy = 400;
    else
        subtitle = 'Percentage of regression (%)';
        ylim([0 70])
        yyy = 57;
    end
    set(gca,'box','off','LineWidth',1)
    set(gca,'XTickLabel',{'Pre-target','Target'},'FontSize',7,'FontWeight','normal','FontName','Arial')
    text(1,yyy,['p = ' num2str(round(1000*tmp.Pre_stat.p)/1000)],'FontSize',7,'FontWeight','normal','FontName','Arial');
    text(2,yyy,['p = ' num2str(round(1000*tmp.Tag_stat.p)/1000)],'FontSize',7,'FontWeight','normal','FontName','Arial');
    title(subtitle,'FontSize',7,'FontWeight','normal','FontName','Arial');
    if ff == 1
        legendflex(a,[{'Incong'},{'Cong'}],'anchor', {'ne','ne'}, 'buffer', [5 -5],'fontsize',7,'FontWeight','normal','FontName','Arial','xscale',0.6,'box', 'off');
    end
    set(gcf, 'renderer', 'painters')
    saveas(gcf,[figpath FigNam{ff}],'svg');
end


%% %%======= violin plot of the word-freq effect on coherence
colmat = [0 114 189;217 83 25]./255;
figtitle = 'Ttest_Freq_violin';
nsigsub = size(TagCoh.PreTarg_Ttest.data4test,1);
group = [cellstr(repmat('Incong',nsigsub,1)); cellstr(repmat('Cong',nsigsub,1))];
grouporder={'Incong','Cong'};
EpochType = {'PreTarg','Targ'};
subtitles = {'Pre-target interval','Target interval'};
figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 450 250]);
for mmm = 1:length(EpochType)
    eval(['vdata = [TagCoh.' EpochType{mmm} '_Ttest.data4test(:,1); TagCoh.' EpochType{mmm} '_Ttest.data4test(:,2)];']);
    h = subplot(1,length(EpochType),mmm);
    vp = violinplot(vdata, group,'GroupOrder',grouporder);
    vp(1).ViolinColor = colmat(1,:);
    vp(2).ViolinColor = colmat(2,:);
    vp(1).ShowMean = 1; vp(2).ShowMean = 1;
    vp(1,1).MedianPlot.Visible = 'off';
    vp(1,2).MedianPlot.Visible = 'off';
    vp(1,1).MeanPlot.LineWidth = 1.5;
    vp(1,2).MeanPlot.LineWidth = 1.5;
    vp(1,1).BoxWidth = 0.01; vp(1,2).BoxWidth = 0.01;
    ylabel('Coherence at 60 Hz (r^2)','FontSize',7,'FontWeight','normal','FontName','Arial');
    xlabel('Semantic congruency of target','FontSize',7,'FontWeight','normal');
    set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
    set(gca,'box','off','LineWidth',1)
    %%% plot the line linking each subject
    hold on;
    x1 = vp(1,1).ScatterPlot.XData;
    y1 = vp(1,1).ScatterPlot.YData;
    x2 = vp(1,2).ScatterPlot.XData;
    y2 = vp(1,2).ScatterPlot.YData;
    plot([x1; x2],[y1; y2],'Color',[.8 .8 .8],'linewidth',0.5)
    %%% add stat
    if mmm == 1
        plot([1 2],[0.065 0.065],'k','LineWidth',1)
        text(1.5,0.061,'*','FontWeight','normal','FontSize',14,'FontName','Arial')
    else
        plot([1 2],[0.065 0.065],'k','LineWidth',1)
        text(1.5,0.063,'n.s.','FontWeight','normal','FontSize',7,'FontName','Arial')
    end
    title(subtitles{mmm},'FontSize',7,'FontWeight','normal','FontName','Arial');
end
set(gcf, 'renderer', 'painters')
saveas(h,figtitle);
saveas(h,[figpath figtitle],'svg');


%% raincloud plot for coherence diff
colmat = [128 128 128]./255;
figtitle = 'Ttest_diff_violin';
figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 170 250]);
diff_pre = TagCoh.PreTarg_Ttest.data4test(:,1)-TagCoh.PreTarg_Ttest.data4test(:,2);
diff_tar = TagCoh.Targ_Ttest.data4test(:,1)-TagCoh.Targ_Ttest.data4test(:,2);
data{1,1} = diff_pre;
data{2,1} = diff_tar;
h = rm_raincloud(data, colmat);
set(gca, 'XLim', [-0.01 0.08]);


%% scatter plot for coherence diff
figtitle = 'Ttest_diff_scatter';
figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 170 250]);
diff_pre = TagCoh.PreTarg_Ttest.data4test(:,1)-TagCoh.PreTarg_Ttest.data4test(:,2);
diff_tar = TagCoh.Targ_Ttest.data4test(:,1)-TagCoh.Targ_Ttest.data4test(:,2);
ydata = [diff_pre diff_tar];
[r, c] = size(ydata);
xdata = repmat(1:c, r, 1);
% for explanation see 
% http://undocumentedmatlab.com/blog/undocumented-scatter-plot-jitter
scatter(xdata(:), ydata(:), 20,[0.5 0.5 0.5], 'jitter','on', 'jitterAmount', 0.05);
hold on;
plot([xdata(1,:)-0.3; xdata(1,:) + 0.3], repmat(mean(ydata, 1), 2, 1), 'k-','LineWidth',1)
ylabel('coherence at 60 Hz (r2, incong-cong)','FontSize',7,'FontWeight','normal','FontName','Arial');
xticks([1 2]);
xticklabels({'Pre-target','target'});
set(gcf, 'renderer', 'painters')
saveas(gcf,figtitle);
saveas(gcf,[figpath figtitle],'svg');


%% bar with error plot for the coherence diff
figtitle = 'Ttest_diff_bar';
figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 170 250]);
diff_pre = TagCoh.PreTarg_Ttest.data4test(:,1)-TagCoh.PreTarg_Ttest.data4test(:,2);
diff_tar = TagCoh.Targ_Ttest.data4test(:,1)-TagCoh.Targ_Ttest.data4test(:,2);
y = [mean(diff_pre); mean(diff_tar);];
errY = [std(diff_pre); std(diff_tar);]./sqrt(size(diff_pre,1));
h = barwitherr(cat(3,-errY,errY), y);    % Plot with errorbars
set(h,'FaceColor',[.5 .5 .5])
set(gca,'box','off','LineWidth',1)
set(gca,'XTickLabel',{'Pre-target','Target'},'FontSize',7,'FontWeight','normal','FontName','Arial')
ylabel('Coherence at 60 Hz (r^2, incong-cong)','FontSize',7,'FontWeight','normal','FontName','Arial');
set(gcf,'renderer', 'painters')
saveas(gcf,figtitle);
saveas(gcf,[figpath figtitle],'svg');


%% === plot the group level curve with sig, only for HM (half-maximum)
%%% PreTarg_HM_tim_tvalue=2.17,threshold for p=0.05 is
%%% t=2.045(two-tail,n=29)
curve_xmax = 0.4; %%% the x-axis range
etp = find(TagCoh.time < curve_xmax); %%% rt window aligned with zero--saccadeonset
zero_tp = nearest(TagCoh.time,0);
tprange = zero_tp:etp(end);
timrange = TagCoh.time(tprange);
figtitles = {'Pre-target','Target'};
EpochType = {'PreTarg','Targ'};
for mmm = 1:length(EpochType)
    eval(['cohdata_all = TagCoh.' EpochType{mmm} '_Coh_TimSubCond(tprange,:,[2 3]);']); %% time*sub*cond
    meanRFT = squeeze(mean(cohdata_all,2));
    seRFT = squeeze(nanstd(cohdata_all,0,2)./sqrt(length(TagCoh.SigSubID)));
    eval(['HM = TagCoh.Jackknife_latency.' EpochType{mmm} '_HM_tim(end,:,1);']);
    colmat = [0 114 189;217 83 25]./255;
    %%%% figure
    figtitle = [figtitles{mmm} ' coherence onset latency'];
    h = figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 220 240]);
    a = shadedErrorBar(timrange,meanRFT(:,1),seRFT(:,1),{'color',colmat(1,:)},0.8);hold on;
    b = shadedErrorBar(timrange,meanRFT(:,2),seRFT(:,2),{'color',colmat(2,:)},0.9);
    a.mainLine.LineWidth = 1;
    b.mainLine.LineWidth = 1;
    legendflex([a.mainLine,b.mainLine],{'Incong';'Cong'},'anchor', {'ne','ne'}, 'buffer', [0 0],'Fontsize',7,'xscale',1,'box','off');
    set(gca,'FontSize',7,'FontWeight','normal');
    set(gca,'box','off','LineWidth',1)
    set(gca,'XTick',0:0.05:0.4);
    set(gca,'XTickLabel',{'FixOn','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4'},'FontSize',7)
    set(gca,'YTick',0.01:0.01:0.06);
    title(figtitle,'FontSize',7,'FontWeight','normal')
    xlabel('Time (s)','FontSize',7,'FontWeight','normal')
    ylabel('Coherence at 60 Hz (r^2)','FontSize',7,'FontWeight','normal','FontName','Arial');
    %%% plot the latency lines
    hold on;
    timeliney = get(gca,'ylim');
    timelinex = repmat(HM(1),1,2);
    plot(timelinex', timeliney','--','color',colmat(1,:),'LineWidth',1);
    timelinex = repmat(HM(2),1,2);
    plot(timelinex', timeliney','--','color',colmat(2,:),'LineWidth',1);
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'svg');
end




%% fixation-related fields -- using sensors from Wang et.al., 2012 or Pan et.al.,2021
cd Z:\Semantic\Results\sv\ERF_FirstFixations\
load('ERF_rmbl_cmb.mat')
sentype = {'Wang2012','Pan2021'};
EpochType = {'PreTarg';'Targ'};
figure('color', [1 1 1],'name','Avg_FRFs','Position',[100 100 450 320]);
for sen = 1:length(sentype)
    load([sentype{sen} '_SigSen.mat']) % Sens
    %%% topo for the sensors
    figtitle = ['sensors from ' sentype{sen}];
    % get a random erf
    erf = ERF.PreTarg_incong{1,1};
    tmp = zeros(size(erf.label));
    for bbb = 1:length(Sens)
        idx = find(strcmp(erf.label,Sens{bbb}));
        tmp(idx) = tmp(idx)+1;
    end
    erf.avg = repmat(tmp,1,size(erf.avg,2));
    cfg =[];
    cfg.layout = 'neuromag306cmb.lay';
    cfg.comment = ' ';
    h = subplot(2,3,1+3*(sen-1));
    ft_topoplotER(cfg,erf);
    title(['Topo for ' figtitle],'FontWeight','normal','FontSize',7);
    colorbar;
    cmap = colormap(flipud(cmap));
    
    %%% plot the FRPs curves over participants
    n_tim = length(erf.time);
    n_sub = length(ERF.subs);
    erf_all = nan(2,n_tim,n_sub);
    for mmm = 1:length(EpochType)
        for s = 1:n_sub
            eval(['tmp_1 = ERF.' EpochType{mmm} '_incong{1,s};']);
            eval(['tmp_2 = ERF.' EpochType{mmm} '_cong{1,s};']);
            %%% remove baseline again after combining planar sensors to get a
            %%% near-zero baseline for the plot
            cfg          = [];
            cfg.baseline = [-0.2 0];
            tmp_1        = ft_timelockbaseline(cfg, tmp_1);
            tmp_2        = ft_timelockbaseline(cfg, tmp_2);
            eval(['label_all = ERF.' EpochType{mmm} '_incong{1,s}.label;']);
            sigchan_id = cellfun(@(x) find(strcmp(x,label_all)),Sens,'Uni',true);
            tmp_1 = mean(tmp_1.avg(sigchan_id,:),1);
            tmp_2 = mean(tmp_2.avg(sigchan_id,:),1);
            erf_all(:,:,s) = [tmp_1; tmp_2];
        end
        % to transfer the unit from T to fT
        erf_all = erf_all.*10e12;
        % plot
        colmat = [0 114 189;217 83 25]./255;
        figtitle = [EpochType{mmm} 'et averaged FRFs'];
        h = subplot(2,3,mmm+1+3*(sen-1));
        meanERF = mean(erf_all,3);
        seERF = nanstd(erf_all,0,3)./sqrt(n_sub);
        a = shadedErrorBar(erf.time,meanERF(1,:),seERF(1,:),{'color',colmat(1,:)},0.8);hold on;
        b = shadedErrorBar(erf.time,meanERF(2,:),seERF(2,:),{'color',colmat(2,:)},0.9);
        a.mainLine.LineWidth = 1;
        b.mainLine.LineWidth = 1;
        set(gca,'box','off','LineWidth',1)
        xlim([-0.2 0.5]);
        set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
        title(figtitle,'FontSize',7,'FontName','Arial','FontWeight','normal');
        xlabel('Time (s)','FontWeight','normal','FontSize',7,'FontName','Arial','FontWeight','normal');
        ylabel('Fixation Related Fields (fT/cm)','FontSize',7,'FontWeight','normal');
    end
    legendflex([a.mainLine,b.mainLine],{'Incong';'Cong'},'anchor', {'ne','ne'}, 'buffer', [0 0],'Fontsize',7,'xscale',1,'box','off');
end
set(gcf, 'renderer', 'painters')
saveas(gcf,'Avg_FRFs');
saveas(gcf,'Avg_FRFs','svg');









%% correlation between Tag/N400 and eye movement matrics
%%%% plotting correlation results using python
% TagPre_N400Targ
x = [Corr.TagPre(:,1)-Corr.TagPre(:,2)]';
y = [Corr.N400Targ(Corr.TagSigSub,1)-Corr.N400Targ(Corr.TagSigSub,2)]';
% copy data into python, place , between every two numbers
sns.jointplot(x, y, kind="reg", truncate=False, color="grey", height = 3)
plt.savefig('U:/writing/Parafoveal_Semantic/figs/TagPre_N400Targ.svg', dpi = 300)
% put coeff and p values
Corr.TagPre_N400Targ_CoefP

% TagPre_DurationPerWrd
y = [Corr.DurationPerwrd(Corr.TagSigSub,1)]';
% put coeff and p values
Corr.TagPre_DurationPerwrd_CoefP


% TagPre_RegdifTarg
y = [Corr.Regres_Targ(Corr.TagSigSub,1)-Corr.Regres_Targ(Corr.TagSigSub,2)]';
% put coeff and p values
Corr.TagPre_RegdifTarg_CoefP


% N400Targ_DurationPerwrd
x = [Corr.N400Targ(:,1)-Corr.N400Targ(:,2)]';
y = [Corr.DurationPerwrd]';
% put coeff and p values
Corr.N400Targ_DurationPerwrd_CoefP


% N400Targ_RegdifTarg
y = [Corr.Regres_Targ(:,1)-Corr.Regres_Targ(:,2)]';
% put coeff and p values
Corr.N400Targ_RegdifTarg_CoefP

























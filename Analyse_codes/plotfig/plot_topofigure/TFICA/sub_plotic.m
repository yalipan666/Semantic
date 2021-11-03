function [] = sub_plotic(hObject,eventdata,filename,ica_type,n_ic,t,f,s_n,P_n,A,t_baseln_lim,channelfile)
% show 1: wave of all trials; 2: TFD; 3: avg wave; 4: topo

%% FFT
L = numel(t);
N_Trials = size(s_n,2);
Fs = 1/(t(2)-t(1));
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Pfft = fft(s_n-repmat(mean(s_n,1),size(s_n,1),1),NFFT,1)/L;
ff = Fs/2*linspace(0,1,NFFT/2);
ff = ff(find(ff>=min(f) & ff<=max(ff)));
Pfft = nanmean(abs(Pfft(find(ff>=min(f) & ff<=max(ff)),:)).^2,2);

figure('unit','normalized','position',[0.1672    0.0688    0.6500    0.8025])
set(gcf,'name',['IC ',num2str(n_ic),' :: ',ica_type,' :: ',filename],'numbertitle','off')

%% PLOT ALL TRIALS
ax_alltrials = subplot(4,3,[2 3]);
imagesc(t,[1:N_Trials],s_n')
axis xy; hold on; box off;
set(gca,'xlim',[min(t),max(t)],'xticklabel',[])
ylabel('Trials')
title(['IC',num2str(n_ic),' Waveforms'],'fontsize',16)
colorbar
set(gca,'position',get(gca,'position')+[0 -0.030 0 0.05])
pos_alltrials = get(gca,'position');
plot([0 0],get(gca,'ylim'),'k','linewidth',2)

%% PLOT AVG ERP
ax_avgtrial = subplot(4,3,[5 6]);
plot(t,nanmean(s_n,2),'linewidth',2)
box off; hold on;
plot([min(t),max(t)],[0 0],'k','linewidth',1,'linestyle','--')
xlabel('Time (sec)'); ylabel('Amplitude');
set(gca,'ylim',[-1 1]*max(abs(get(gca,'ylim'))))
plot([0 0],get(gca,'ylim'),'k','linewidth',2)
set(gca,'xlim',[min(t),max(t)],'ytick',0.5*[get(gca,'ylim')])
pos_avgtrial = get(gca,'position');
set(gca,'position',[pos_avgtrial(1),pos_alltrials(2)-pos_alltrials(4)/1.8,pos_alltrials(3),pos_alltrials(4)/1.8])

%% PLOT TFR
ax_avgtfr = subplot(4,3,[8 9 11 12]);
% load tfd_cmap.mat
% baseline correction
t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2)));
P_Baseline_Mean_Vec = mean(P_n(:,t_baseln_idx,:),2);
P_Baseline_Mean = repmat(P_Baseline_Mean_Vec,[1,size(P_n,2),1]);
P_BC = (P_n - P_Baseline_Mean );
P_BC_AVG = nanmean((P_BC),3);

clim = [-1 1]*max(abs(P_BC_AVG(:)));
imagesc(t,f,P_BC_AVG,clim)
% set(gcf,'colormap',tfd_cmap)
hold on;
axis xy;set(gca,'xlim',[min(t),max(t)])
plot([0 0],get(gca,'ylim'),'k','linewidth',2)
xlabel('Time (sec)')
hc = colorbar;
title('Average Time-frequency Representation','fontsize',16)

pos_avgtfr = get(ax_avgtfr,'position');
gca_pbtfr = [pos_avgtfr(1)+pos_avgtfr(3)*1.08,pos_avgtfr(2)+pos_avgtfr(4)*1.1,pos_avgtfr(3)/3.2,pos_avgtfr(4)/5];
pb_alltfr = uicontrol('Style', 'pushbutton', 'String',['Show All TFRs'],'units','normalized',...
        'fontsize',10, 'BackgroundColor','y','ForegroundColor','r','FontWeight','bold',...
        'Position', gca_pbtfr, 'Callback', {'sub_plotalltfr',filename,ica_type,n_ic,t,f,P_n,t_baseln_lim});

gca_fig1 = get(hc,'position');
pb_tfrpower = uicontrol('Style', 'pushbutton', 'String',['Power'],'units','normalized','BackgroundColor','b','ForegroundColor','w','FontWeight','bold',...
        'Position', gca_fig1+[gca_fig1(3)*2 gca_fig1(4)*0.80 gca_fig1(3) -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'power','avgtrials'});
pb_tfrdiff = uicontrol('Style', 'pushbutton', 'String',['DIFF'],'units','normalized','BackgroundColor','r','ForegroundColor','w','FontWeight','bold',...
        'Position', gca_fig1+[gca_fig1(3)*2 gca_fig1(4)*0.55 gca_fig1(3) -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'diff','avgtrials'});
pb_tfrer = uicontrol('Style', 'pushbutton', 'String',['ER%'],'units','normalized','BackgroundColor','g','ForegroundColor','w','FontWeight','bold',...
        'Position', gca_fig1+[gca_fig1(3)*2 gca_fig1(4)*0.30 gca_fig1(3) -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'er','avgtrials'});
pb_tfrlogr = uicontrol('Style', 'pushbutton', 'String',['Log%'],'units','normalized','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',...
        'Position', gca_fig1+[gca_fig1(3)*2 gca_fig1(4)*0.05 gca_fig1(3) -gca_fig1(4)*0.75], 'fontsize',8, 'horizontalalignment','right','Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'log','avgtrials'});


%% PLOT SPECTRUM
subplot(4,3,[7 10])
% plot(ff,10*log10(Pfft),'k','linewidth',2)
plot(ff,(Pfft),'k','linewidth',2)
set(gca,'View',[-90 90],'xlim',[min(f) max(f)])
xlabel('Frequency (Hz)')
% ylabel('Power (10*log_1_0(\it\mu\rmV^2/Hz))')
ylabel('Power (\it\mu\rmV^2/Hz)')
title('Average Spectrum','fontsize',16)

%% PLOT TOPOGRAPHY
subplot(4,3,[1 4])
eeglab_topoplot(A(:,n_ic),channelfile);   
title(['IC',num2str(n_ic),' Topography'],'fontsize',16)
function [] = sub_plotalltfr( hObject,eventdata,filename,ica_type,n_ic,t,f,P_n, t_baseln_lim )
% plot all TFRs

%% determine the row and column in a figure
N_Trials = size(P_n,3);
n_row = ceil(sqrt(N_Trials));
n_col = ceil(sqrt(N_Trials));
if (n_row-1)*n_col>=N_Trials; n_row=n_row-1; end;

t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2)));
P_Baseline_Mean_Vec = mean(P_n(:,t_baseln_idx,:),2);
P_Baseline_Mean = repmat(P_Baseline_Mean_Vec,[1,size(P_n,2),1]);
P_BC = (P_n - P_Baseline_Mean );

%%
figure('unit','normalized','position',[0    0.0375    1.0000    0.8700],'Toolbar','figure')
set(gcf,'name',['TFRs: IC ',num2str(n_ic),' :: ',ica_type,' :: ',filename],'numbertitle','off')
% load tfd_cmap.mat
for n_trial=1:N_Trials
    ax = subplot(n_row,n_col,n_trial);
    P_show = P_BC(:,:,n_trial);
    clim = [-1 1]*max(abs(P_show(:)));
    imagesc(t,f,P_show,clim)
%     set(gcf,'colormap',tfd_cmap)
    axis xy; hold on;
    set(gca,'fontsize',8)
    plot([0 0],get(gca,'ylim'),'k','linewidth',2)
    title(['Trial ',num2str(n_trial)],'fontsize',get(gca,'fontsize')+2)
    if mod(n_trial,n_col)==1
        ylabel('Frequency (Hz)')
    end
    if n_trial>(n_row-1)*n_col
        xlabel('Time (sec)')
    end
    if n_trial==N_Trials
        hc = colorbar('east');
        pos_hc = get(hc,'position');
        pos_ax = get(ax,'position');
        set(hc,'position',[pos_ax(1)+pos_ax(3)*1.1 pos_ax(2) pos_hc(3) pos_ax(4)])
    end
    if n_trial==1
        gca_fig1 = get(gca,'position');
        pb_tfrpower = uicontrol('Style', 'pushbutton', 'String',['Power'],'units','normalized','BackgroundColor','b','ForegroundColor','w','FontWeight','bold',...
            'Position', gca_fig1+[-gca_fig1(3) gca_fig1(4)*0.80 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'power','alltrials'});
        pb_tfrdiff = uicontrol('Style', 'pushbutton', 'String',['DIFF'],'units','normalized','BackgroundColor','r','ForegroundColor','w','FontWeight','bold',......
            'Position', gca_fig1+[-gca_fig1(3) gca_fig1(4)*0.55 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'diff','alltrials'});
        pb_tfrer = uicontrol('Style', 'pushbutton', 'String',['ER%'],'units','normalized','BackgroundColor','g','ForegroundColor','w','FontWeight','bold',......
            'Position', gca_fig1+[-gca_fig1(3) gca_fig1(4)*0.30 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'er','alltrials'});
        pb_tfrlogr = uicontrol('Style', 'pushbutton', 'String',['Log%'],'units','normalized','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',......
            'Position', gca_fig1+[-gca_fig1(3) gca_fig1(4)*0.05 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'fontsize',8, 'horizontalalignment','right','Callback', {'sub_plottfr',P_n,t,f,t_baseln_lim,'log','alltrials'});
    end
end

function sub_plottfr(hObject,eventdata,P,t,f,t_baseln_lim,tfr_opt,func_opt)
% to show TFR in different manner

NDims = ndims(P); % dimension of P

if NDims==4 % chan*freq*time*trial %plot ic or trial
    N_Comp = size(P,1); % or channel
    % determine the row and column in a figure
    n_row = ceil(sqrt(N_Comp));
    n_col = ceil(sqrt(N_Comp));
    if (n_row-1)*n_col>=N_Comp; n_row=n_row-1; end;

    %% Baseline correction
    t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2)));
    % baseline value is obtained from each trial, NOT all trials
    P_Baseline_Mean_Vec = mean(P(:,:,t_baseln_idx,:),3);
    P_Baseline_Mean = repmat(P_Baseline_Mean_Vec,[1,1,size(P,3),1]);
    P_DIFF = (P - P_Baseline_Mean );
    P_ER = P_DIFF./P_Baseline_Mean;

    for n_ic=1:size(P,1)
        
        subplot(n_row,n_col,n_ic);
        if strcmp(tfr_opt,'power')
            P_show = nanmean(squeeze(P(n_ic,:,:,:)),3);            
        elseif strcmp(tfr_opt,'diff')
            P_show = nanmean(squeeze(P_DIFF(n_ic,:,:,:)),3);
        elseif strcmp(tfr_opt,'er')
            P_show = nanmean(squeeze(P_ER(n_ic,:,:,:)),3);
        elseif strcmp(tfr_opt,'log')
            P_show = 20*log10(nanmean(squeeze(P_ER(n_ic,:,:,:)),3)+1);        
        end

        clim = [-1 1]*max(abs(P_show(:)));
        imagesc(t,f,P_show,clim)
%         set(gcf,'colormap',tfd_cmap)
        axis xy; hold on
        plot([0 0],get(gca,'ylim'),'k','linewidth',2)
        if mod(n_ic,n_col)==1
            ylabel('Freq (Hz)')
        end
        if n_ic>(n_row-1)*n_col
            xlabel('Time (sec)')
        end
        axis tight
        hold off
        drawnow
    end
elseif (NDims==3)&(strcmp(func_opt,'avgtrials')) % freq*time*trial %plot avg of one ic across trials
    ax_avgtfr = subplot(4,3,[8 9 11 12]);    
    
    % baseline correction
    t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2)));
    P_Baseline_Mean_Vec = mean(P(:,t_baseln_idx,:),2);
    P_Baseline_Mean = repmat(P_Baseline_Mean_Vec,[1,size(P,2),1]);
    P_DIFF = (P - P_Baseline_Mean );
    P_ER = P_DIFF./P_Baseline_Mean;
    if strcmp(tfr_opt,'power')
        P_show = nanmean(squeeze(P),3);
    elseif strcmp(tfr_opt,'diff')
        P_show = nanmean(squeeze(P_DIFF),3);
    elseif strcmp(tfr_opt,'er')
        P_show = nanmean(squeeze(P_ER),3);
    elseif strcmp(tfr_opt,'log')
        P_show = 20*log10(nanmean(squeeze(P_ER),3)+1);
    end

    clim = [-1 1]*max(abs(P_show(:)));
    imagesc(t,f,P_show,clim)
    % set(gcf,'colormap',tfd_cmap)
    hold on;
    axis xy;set(gca,'xlim',[min(t),max(t)])
    plot([0 0],get(gca,'ylim'),'k','linewidth',2)
    xlabel('Time (sec)')
    hc = colorbar;
    title('Average Time-frequency Representation','fontsize',16)

elseif (NDims==3)&(strcmp(func_opt,'alltrials')) % freq*time*trial %plot all trials of one ic
    
    N_Trials = size(P,3);
    n_row = ceil(sqrt(N_Trials));
    n_col = ceil(sqrt(N_Trials));
    if (n_row-1)*n_col>=N_Trials; n_row=n_row-1; end;

    t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2)));
    P_Baseline_Mean_Vec = mean(P(:,t_baseln_idx,:),2);
    P_Baseline_Mean = repmat(P_Baseline_Mean_Vec,[1,size(P,2),1]);
    P_DIFF = (P - P_Baseline_Mean );
    P_ER = P_DIFF./P_Baseline_Mean;

    for n_trial=1:N_Trials
        ax = subplot(n_row,n_col,n_trial);
        if strcmp(tfr_opt,'power')
            P_show = P(:,:,n_trial);
        elseif strcmp(tfr_opt,'diff')
            P_show = P_DIFF(:,:,n_trial);
        elseif strcmp(tfr_opt,'er')
            P_show = P_ER(:,:,n_trial);            
        elseif strcmp(tfr_opt,'log')
            P_show = 20*log10(P_ER(:,:,n_trial)+1);
        end
        
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
%         if n_trial==N_Trials
%             hc = colorbar('east');
%             pos_hc = get(hc,'position');
%             pos_ax = get(ax,'position');
%             set(hc,'position',[pos_ax(1)+pos_ax(3)*1.1 pos_ax(2) pos_hc(3) pos_ax(4)])
%         end
    end
end % end of IF

end % end of function

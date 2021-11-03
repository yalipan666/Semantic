function [] = sub_plotallic(hObject,eventdata,filename,ica_type,var_type,channelfile)
% show eeg/tfr/topo of all ICs

%% variables
global TF_IC T_IC PTF_IC PT_IC t f ChannelLabels t_baseln_lim pathname filename
warning off last
if strcmp(ica_type,'Fourier-ICA')
    s = TF_IC.s;
    A = TF_IC.A;
    A=real(A);
    P = TF_IC.P;
elseif strcmp(ica_type,'Infomax-ICA')
    s = T_IC.s;
    A = T_IC.A;
    P = T_IC.P;
elseif strcmp(ica_type,'Probabilistic Fourier-ICA')
    s = PTF_IC.s;
    A = PTF_IC.A;
    P = PTF_IC.P;
elseif strcmp(ica_type,'Probabilistic Infomax-ICA')
    s = PT_IC.s;
    A = PT_IC.A;
    P = PT_IC.P;
end

N_Trials = size(s,3);
N_Comp = size(P,1); % <==== if #components is diffrent with #channels, it should be modified

%% figure maps
% determine the row and column in a figure
n_row = ceil(sqrt(N_Comp));
n_col = ceil(sqrt(N_Comp));
if (n_row-1)*n_col>=N_Comp; n_row=n_row-1; end;

figure('unit','normalized','position',[0    0.0375    1.0000    0.8700],'Toolbar','figure')
set(gcf,'name',[upper(var_type),' :: ',ica_type,' :: ',filename],'numbertitle','off')
for n_ic=1:N_Comp
    sub_fig(n_ic) = subplot(n_row,n_col,n_ic);
	P_n = squeeze(P(n_ic,:,:,:));
    s_n = squeeze(s(n_ic,:,:));
	if strcmp(var_type,'topo')
        eeglab_topoplot(A(:,n_ic),channelfile);
        axis tight
    elseif strcmp(var_type,'tfr')
%         load tfd_cmap.mat
        
        t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2)));
        P_Baseline_Mean_Vec = mean(P_n(:,t_baseln_idx,:),2);
        P_Baseline_Mean = repmat(P_Baseline_Mean_Vec,[1,size(P_n,2),1]);
        P_BC = (P_n - P_Baseline_Mean);
        
        P_show = nanmean((P_BC),3);
        clim = [-1 1]*max(abs(P_show(:)));
        imagesc(t,f,P_show,clim)
        axis xy; hold on
%         set(gcf,'colormap',tfd_cmap)
        plot([0 0],get(gca,'ylim'),'k','linewidth',2)
        if mod(n_ic,n_col)==1
            ylabel('Freq (Hz)')
        end
        if n_ic>(n_row-1)*n_col
            xlabel('Time (sec)')
        end
        axis tight
        if n_ic==1
            gca_fig1 = get(sub_fig(n_ic),'position');
            pb_tfrpower = uicontrol('Style', 'pushbutton', 'String',['Power'],'units','normalized','BackgroundColor','b','ForegroundColor','w','FontWeight','bold',...
            'Position', gca_fig1+[-gca_fig1(3)*.8 gca_fig1(4)*0.8 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P,t,f,t_baseln_lim,'power'});
            pb_tfrdiff = uicontrol('Style', 'pushbutton', 'String',['Diff'],'units','normalized','BackgroundColor','r','ForegroundColor','w','FontWeight','bold',......
            'Position', gca_fig1+[-gca_fig1(3)*.8 gca_fig1(4)*0.55 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P,t,f,t_baseln_lim,'diff'});
            pb_tfrer = uicontrol('Style', 'pushbutton', 'String',['ER%'],'units','normalized','BackgroundColor','g','ForegroundColor','w','FontWeight','bold',......
            'Position', gca_fig1+[-gca_fig1(3)*.8 gca_fig1(4)*0.3 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'horizontalalignment','right','fontsize',8, 'Callback', {'sub_plottfr',P,t,f,t_baseln_lim,'er'});
            pb_tfrlogr = uicontrol('Style', 'pushbutton', 'String',['Log%'],'units','normalized','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',......
            'Position', gca_fig1+[-gca_fig1(3)*.8 gca_fig1(4)*0.05 -gca_fig1(3)/2 -gca_fig1(4)*0.75], 'fontsize',8, 'horizontalalignment','right','Callback', {'sub_plottfr',P,t,f,t_baseln_lim,'log'});
        end           
    elseif strcmp(var_type,'eeg')
        scale_eeg = round(N_Trials/4)/(max(max(abs(mean(s(:,:,:),3))))); %universal scale for channels
        %scale_eeg = round(N_Trials/4)/(max(abs(mean(s(n_ic,:,:),3)))); %individual scale based on average
        %scale_eeg = round(N_Trials/4)/(max(max(abs(s(n_ic,:,:))))); %individual scale based on all
        imagesc(t,[1:N_Trials],s_n',[-1 1]*(max(max(abs(s(n_ic,:,:))))))
        axis xy; hold on;
        ytick = get(gca,'ytick');
        ytick = ytick(find(ytick>=0));
        set(gca,'ytick',ytick);
        set(gca,'ylim',[-round(N_Trials/2),N_Trials])
        plot(t,nanmean(s_n,2)*scale_eeg-round(N_Trials/4),'b','linewidth',1)
        plot([0 0],get(gca,'ylim'),'k','linewidth',2)
        if mod(n_ic,n_col)==1
            ylabel('Trial')
        end
        if n_ic>(n_row-1)*n_col
            xlabel('Time (sec)')
        end
    end

    set(gca,'fontsize',6)
	gca_pos = get(gca,'position');
    gca_x = gca_pos(1);    gca_y = gca_pos(2);
    gca_l = gca_pos(3);    gca_h = gca_pos(4); 
    gca_pb(1) = gca_x+gca_l/4;
    gca_pb(2) = gca_y+gca_h;
    gca_pb(3) = gca_l/2;
    gca_pb(4) = gca_h/4;
    pb = uicontrol('Style', 'pushbutton', 'String',['IC',num2str(n_ic)],'units','normalized',...
        'Position', gca_pb, 'fontsize',8, 'Callback', {'sub_plotic',filename,ica_type,n_ic,t,f,s_n,P_n,A,t_baseln_lim,channelfile});
    drawnow
end
end % end of function


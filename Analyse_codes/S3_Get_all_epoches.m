% % % copy from Lexical/Analyse_codes
% % % 20210719 remove artefcts based on ica components, get all kinds of epochs based on
% % % event

function S3_Get_all_epoches(server,sid)
TaskId = 1; %[1 3]; DataSets = {'sv','of','fa'}; 'sv and of' tasks are combined in 1st task

%%% set paths
if server
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
    rootdir = '/rds/projects/2018/jenseno-reading/Semantic/';
else
    addpath Z:\fieldtrip-20200220\
    ft_defaults
    rootdir = 'Z:\Semantic\';
end

%%% basic settingup
RunCond = 'WrdOn'; %%'WrdOn';%%% epoch aligned with RunCond
PPara.filename = RunCond;
AnaTW = 1000;
epoch_select = 'abs(PPara.event_all(:,3))<2'; %% select: n-1,n,n+1
DoBaseline = 0; % get epochs for baseline period
ica = 0; % remove ica components
%%% get file names
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
DataSets = ExpInfo.DSName;

for ddd = TaskId
    
    DS = DataSets{ddd};
    conds = ExpInfo.CondID;% condition mat for diff tasks
    eval(['subjects = ExpInfo.subjects.' DS ';']);
    sub = subjects{sid};
    PPara.sub = sub;
    PPara.badsens = ExpInfo.BadSensor{sid};
    PPath.SaveData = [rootdir 'Analyse_data' filesep DS '_' sub filesep];
    
    %% ========== artefacts removal based on ICA components after preprocessing
    if ica
        load([PPath.SaveData 'ica']); %'data4ICA','comp','Trig'
        
        %%% plot the components for visual inspection
        figure
        cfg           = [];
        cfg.component = 1:30;       % specify the component(s) that should be plotted
        cfg.layout    = 'neuromag306mag.lay'; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp)
        colormap jet;
        
        %%% identifying components
        cfg            = [];
        cfg.channel    = 1:15;
        cfg.continuous = 'no';
        cfg.viewmode   = 'component';
        cfg.layout     = 'neuromag306mag.lay';
        ft_databrowser(cfg, comp);
        colormap jet;
        
        %%% rejecting components back to the un-down-sampling data
        cfg = [];
        cfg.component = input('component(s) ID that needed to be removed []: ');
        data4ICA = ft_rejectcomponent(cfg, comp, data4ICA);
        clear comp
        
        %%% put the ica clean epoch data back to the raw data
        load([PPath.SaveData 'data']); %raw meg data
        load([PPath.SaveData 'hdr']); %meg hdr
        validmegid = cellfun(@(x) find(strcmp({x},hdr.label)),data4ICA.label);
        for i = 1:length(Trig.fixon)
            data(validmegid,Trig.fixon(i):Trig.sentoff(i)) = data4ICA.trial{i};
        end
        save([PPath.SaveData 'data_icaclean.mat'], 'data','-v7.3')
        clear data4ICA
        
        % delete ica.mat
        delete([PPath.SaveData 'ica.mat']); %'data4ICA','comp','Trig'
    else
        load([PPath.SaveData 'data_icaclean.mat']); % data
        load([PPath.SaveData 'hdr']); %meg hdr
    end
    
    %% ================epoching
    load([PPath.SaveData 'Event'])
    %%% get trialinfo index
    PPara.SR = 1000; %sampling rate
    PPara.timezero = 'fixation_on_MEG';
    trig_col = find(strcmp(Event.event_raw_header,PPara.timezero));
    fixdur_col = find(strcmp(Event.event_raw_header,'fixation_duration'));
    cond_col = find(strcmp(Event.event_raw_header,'SentenceCondition'));
    
    %%%% get the baseline epochs
    if DoBaseline
        TrigOn = Event.Trigger_MEG(Event.Trigger_MEG(:,1) == ExpInfo.Trigger.Fix,2);
        epoch_leng = AnaTW;
        event_all = nan(length(TrigOn),size(Event.event_raw,2));
        event_all(:,trig_col) = TrigOn;
        event_all(:,fixdur_col) = epoch_leng.*ones(size(event_all,1),1);
        PPara.event_all = event_all;
        PPara.pretrig = 0;
        PPara.posttrig = PPara.pretrig + epoch_leng; %epoch end timepoint, aligned with marker
        epoch_BL_Cross = Get_Epoch(hdr,data,PPara,trig_col);
        save([PPath.SaveData 'epoch_BL_Cross'],'epoch_BL_Cross','-v7.3');
    end
    
    %%% get the epochs for sentence reading
    PPara.pretrig = -AnaTW/2;
    PPara.posttrig = AnaTW/2;
    tw = [80 AnaTW];
    validduration = Event.event_raw(:,fixdur_col)>= tw(1) & Event.event_raw(:,fixdur_col)<= tw(2);
    PPara.event_all = Event.event_raw(validduration,:);
    eval(['PPara.event_all = PPara.event_all(' epoch_select ',:);']);
    %%%% change the fixation onset trigger to fixation offset
    if contains(RunCond,'Off')
        PPara.event_all(:,trig_col) = PPara.event_all(:,trig_col) + PPara.event_all(:,fixdur_col);
    end
    epoch = Get_Epoch(hdr,data,PPara,trig_col);
    
    %%% seperate epochs into 2 tasks and save them out
    if ddd == 3
        save([PPath.SaveData 'epoch_' PPara.filename],'epoch','-v7.3');
    else %%% seperate epochs into 2 tasks and save them out
        epoch_all = epoch;
        % first task
        tmptrl = [];
        tmpc = conds{1};
        for ccc = 1:length(tmpc)
            tmptrl = [tmptrl epoch_all.trialinfo(:,cond_col)==tmpc(ccc)];
        end
        tmptrl = sum(tmptrl,2);
        cfg = [];
        cfg.trials = find(tmptrl~=0);
        epoch = ft_selectdata(cfg,epoch_all);
        save([PPath.SaveData 'epoch_' PPara.filename],'epoch','-v7.3');
        
        % second task
        tmptrl = [];
        tmpc = conds{2};
        for ccc = 1:length(tmpc)
            tmptrl = [tmptrl epoch_all.trialinfo(:,cond_col)==tmpc(ccc)];
        end
        tmptrl = sum(tmptrl,2);
        cfg = [];
        cfg.trials = find(tmptrl~=0);
        epoch = ft_selectdata(cfg,epoch_all);
        path_task_2 = [rootdir 'Analyse_data' filesep DataSets{2} '_' sub filesep];
        mkdir(path_task_2);
        save([path_task_2 'epoch_' PPara.filename],'epoch','-v7.3');
        if DoBaseline
            save([path_task_2 'epoch_BL_Cross'],'epoch_BL_Cross','-v7.3');
        end
    end
    save([rootdir 'Analyse_data' filesep 'ExpInfo.mat'],'ExpInfo');
    disp(['*** epoching done! ' DS '---' sub]);
    clear data
end




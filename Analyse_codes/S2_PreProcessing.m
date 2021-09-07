% copy from Lexical/Analyse_codes
% 20210720 pre-processing data on bluebear: get raw megdata, eyedata, event,
% and all ica components

function S2_PreProcessing(sid) %id of sub
% number of .fif files in SV task
nfif_sv = 3; %normally 2 .fif files, 3 files for 20210826_b38f

%%% set paths
server = 0;
if server
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
    rootdir = '/rds/projects/2018/jenseno-reading/Semantic/';
else
    addpath Z:\fieldtrip-20200220\
    ft_defaults
    rootdir = 'Z:\Semantic\';
end

PPath.RawMEG = [rootdir 'RawData' filesep 'MEG_data' filesep];
PPath.RawPTB = [rootdir 'RawData' filesep 'PTB_data' filesep];
PPath.RawEye = [rootdir 'RawData' filesep 'EyeLink_data' filesep];

%%% settings for pre-processing
PPara.bpfilter   = 'yes';
PPara.bpfreq     = [0.5 100];
PPara.detrend    = 'yes';
PPara.SR         = 1000;

%%% get file names
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
DataSets = ExpInfo.DSName;

for ddd = [1 3]% DataSets = {'sv','of','fa'}; the data of 'sv and of' tasks are combined in a single dataset 
    % since 'sv, of' tasks are combine in the fist dataset, there's no need
    % to run ddd = 2
    DS = DataSets{ddd};
    eval(['subjects = ExpInfo.subjects.' DS ';']);
    eval(['PTBFiles = ExpInfo.PTBFiles.' DS ';']);
    eval(['EyeFiles = ExpInfo.EyeFiles.' DS ';']);
    sub = subjects{sid};
    PPara.sub = sub;
    PPara.badsens = ExpInfo.BadSensor{sid};
    File.PTB = PTBFiles{sid};
    File.Eye = EyeFiles{sid};
    PPath.SaveData = [rootdir 'Analyse_data' filesep DS '_' sub filesep];
    if ~exist(PPath.SaveData,'dir')
        mkdir(PPath.SaveData);
    end
    
    %%% display processing sub
    disp(['*** PreProcessing: ' DS '--sub--' sub]);
    
    %% get raw MEG data
    if ~exist([PPath.SaveData 'data.mat'],'file')
        if strcmp(DS,'sv') % SV task with two/three .fif files
            tmpf = [sub filesep sub(3:8) filesep sub(10:end)];
            File.MEG = cell(1,nfif_sv);
            for f = 1:nfif_sv
                if f == 1
                   File.MEG{f} = [tmpf '.fif'];
                else
                   File.MEG{f} = [tmpf '-' num2str(f-1) '.fif'];
                end
            end
        else % another task with one .fif file
            File.MEG = {[sub filesep sub(3:8) filesep DS '_' sub(10:end) '.fif']};
        end
        [hdr,data,Trigger_MEG] = Get_MEGData(PPath,File);
        save([PPath.SaveData 'hdr'],'hdr');
        save([PPath.SaveData 'Trigger_MEG'],'Trigger_MEG','-v7.3');
        save([PPath.SaveData 'data'],'data','-v7.3');
    else
        load([PPath.SaveData 'hdr']);
        load([PPath.SaveData 'data']);
        load([PPath.SaveData 'Trigger_MEG']);
    end
    
    %% remove artefacts with ICA on server
    % just run ica for the exp data, excluding the hpi signal at the begining
    % and the end of exp
    PPara.icastart = Trigger_MEG(1,2);
    PPara.cutlen4ica = 2^(nextpow2(4*PPara.SR));%random length for epoching in ica
    PPara.icaend = Trigger_MEG(end,2)+PPara.cutlen4ica;
    [data4ICA, comp] = ICA4rawdata(PPara,hdr,data);
    save([PPath.SaveData 'ica.mat'],'data4ICA','comp','-v7.3')
    
    %% get eyemovement metrics from eyelink
    %first converting data from .edf to .asc (C:/toolbox/SR Research/edfconverter/)
    eyefile = [PPath.RawEye File.Eye '.asc'];
    load([PPath.RawPTB File.PTB],'Result');
    WordLocMat = 2.*Result.WordLocation;
    EyeData = Get_EyeData(eyefile, WordLocMat,ExpInfo.Trigger);
    save([PPath.SaveData 'EyeData'],'EyeData');
    
    %% get event
    load([PPath.RawPTB File.PTB],'Para'); % Para
    Event = Get_Event(EyeData,Para.CondMat,Para.TargLoc,Trigger_MEG,ExpInfo.Trigger,ExpInfo.EventHdr);
    save([PPath.SaveData 'Event'],'Event','-v7.3');
    
    disp(['*** PreProcessing done! ' DS '---' sub]);
end





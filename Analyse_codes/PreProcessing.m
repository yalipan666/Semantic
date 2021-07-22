% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise
% 20210720 pre-processing data on bluebear: get raw megdata, eyedata, event,
% ica components

function PreProcessing(sss)
ddd = 3;

%%% set paths
server = 1;
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
DataSets = {'sv','of','fa'};
DS = DataSets{ddd};
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
eval(['subjects = ExpInfo.subjects.' DS ';']);
eval(['PTBFiles = ExpInfo.PTBFiles.' DS ';']);
eval(['EyeFiles = ExpInfo.EyeFiles.' DS ';']);
sub = subjects{sss};
PPara.sub = sub;
PPara.badsens = ExpInfo.BadSensor{sss};
File.PTB = PTBFiles{sss};
File.Eye = EyeFiles{sss};
PPath.SaveData = [rootdir 'Analyse_data' filesep DS '_' sub filesep];
if ~exist(PPath.SaveData,'dir')
    mkdir(PPath.SaveData);
end

%%% display processing sub
disp(['*** PreProcessing: ' DS '--sub--' sub]);

%% get raw MEG data
if ~exist([PPath.SaveData 'data.mat'],'file')%due to connect unexpected broken data
    if strcmp(DS,'sv') % SV task with two .fif files
        tmpf = [sub filesep sub(3:8) filesep sub(10:end)];
        File.MEG = {[tmpf '.fif'],[tmpf '-1.fif']};
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
CondMat = Para.CondMat;
Event = Get_Event(EyeData,CondMat,Trigger_MEG,ExpInfo.Trigger);
save([PPath.SaveData 'Event'],'Event','-v7.3');

disp(['*** PreProcessing done! ' DS '---' sub]);
end




% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise

function [hdr,data,Trigger_MEG] = Get_MEGData(PPath,File)
%%% read in two MEG datasets
nfif = length(File.MEG);
for p = 1:nfif
    cfg.dataset = [PPath.RawMEG File.MEG{p}];
    eval(['data'  num2str(p) '= ft_read_data(cfg.dataset);']);
    event = ft_read_event(cfg.dataset);
    Trigger_MEG = [[event(strcmp('Trigger',{event.type})).value]' [event(strcmp('Trigger',{event.type})).sample]'];
    eval(['Trigger_MEG' num2str(p) '= Trigger_MEG;']);
end
if nfif == 2 % concatenating 2 datasets
    data_length_1 = size(data1,2);
    Trigger_MEG = [Trigger_MEG1; [Trigger_MEG2(:,1) Trigger_MEG2(:,2)+data_length_1]];
    data = [data1 data2];
    clear data1 data2
else
    data = data1; 
    clear data1
    Trigger_MEG = Trigger_MEG1;
end
hdr =  ft_read_header([PPath.RawMEG File.MEG{1}]);


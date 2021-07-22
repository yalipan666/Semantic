% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise
%store information of the paths, subjects, files, and triggers.

%% subjects and files
DataSets = {'sv','of','fa'};% SV+OF is the semantic task, FA is JEP task
for ddd = 1:length(DataSets)
    if ddd == 1
        % semantic violation sentences
        subjects = {'20210625_b4d4';'20210702_b4bb';'20210702_b5f3';...
            '20210708_b398';'20210708_b5f6';'20210716_b5e6';'20210719_b395';...
            '20210719_b4bc'}; 
        PTBFiles = {'20210625_b4d4';'20210702_b4bb';'20210702_b5f3';...
            '20210708_b398';'20210708_b5f6';'20210716_b5e6';'20210719_b395';...
            '20210719_b4bc'}; 
        EyeFiles = {'0625b4d4';'0702b4bb';'0702b5f3';'0708b398';'0708b5f6';...
            '0716b5e6';'0719b395';'0719b4bc'};
    elseif ddd == 2
        % orth_freq sentences
        subjects = {'20210708_b5f6';'20210716_b5e6';'20210719_b395';...
            '20210719_b4bc'}; 
        PTBFiles = {'of_20210708_b5f6';'20210716_b5e6';'20210719_b395';...
            '20210719_b4bc'}; 
        EyeFiles = {'of_0708b5f6';'0716b5e6';'0719b395';'0719b4bc'};
    else % full_attention_dynamics sentences
        subjects = {'20210625_b4d4';'20210702_b4bb';'20210702_b5f3';...
            '20210708_b398';'20210708_b5f6';'20210716_b5e6';'20210719_b395';...
            '20210719_b4bc'}; 
        PTBFiles = {'20210625_b4d4';'20210702_b4bb';'20210702_b5f3';...
            '20210708_b398';'20210708_b5f6';...
            'fa_20210716_b5e6';'fa_20210719_b395';'fa_20210719_b4bc'}; 
        EyeFiles = {'0625b4d4';'0702b4bb';'0702b5f3';'0708b398';'0708b5f6';...
            'fa_0716b5e6';'fa_0719b395';'fa_0719b4bc'};
    end
    eval(['ExpInfo.subjects.' DataSets{ddd} '= subjects;']);
    eval(['ExpInfo.PTBFiles.' DataSets{ddd} '= PTBFiles;']);
    eval(['ExpInfo.EyeFiles.' DataSets{ddd} '= EyeFiles;']);
end
ExpInfo.BadSensor = {{'1043'};{'1043','1431','2211'};{'1043','2211'};{'1043'};...
    {'1043'};{'1043'};{'1021','1043','2141','2331','2131'};{'1043'}};

% %% ==== gettting names for MRI images: MEGcode_MRIcode
% MEGcode = {'20190531_b588';};
% MRIcode = {'20190606#C467';};
% 
% ExpInfo.MEGcode_MRIcode = {};
% for ccc = 1:length(MEGcode)
%     ExpInfo.MEGcode_MRIcode(ccc,:) = [MEGcode(ccc,1) MRIcode(ccc,1)];
% end

%% triggers
Trigger.Fix = 1;   % fixation onset
Trigger.StartBox = 2; % start box onset
Trigger.SentOn = 4;   %  sentence onset
Trigger.SentOff = 8;  % sentence offset
Trigger.ITI = 16;   % ITI onset
Trigger.OpenEye = 64;   
Trigger.CloseEye = 128;  
Trigger.PureTagOn = 32;   % tagging on in the pure tagging 
ExpInfo.Trigger = Trigger;
save('Z:\Semantic\Analyse_data\ExpInfo.mat','ExpInfo')

















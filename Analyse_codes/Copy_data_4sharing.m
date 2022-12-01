%% copy necessary data for sharing
load('Z:\Semantic\Analyse_data\ExpInfo.mat'); 
! Please delete any unneccessary information (i.e., info about the other two tasks),then save it out to 'Z:\Semantic\DataSharing\'

rawpath = 'Z:\Semantic\RawData\';
dspath = 'Z:\Semantic\DataSharing\';

in_files = {'MEG_data','EyeLink_data','PTB_data'};
out_files = {'MEGdata','Eyelinkdata','PTBdata'};
vars = {'subjects','EyeFiles','PTBFiles'};
ext = {'','.edf','.mat'};
datasets = {'sv'};
nsub = 34; 
for vvv = 1:length(in_files)
    mkdir(dspath,out_files{vvv})
    for ddd = 1:length(datasets)
        for sss = 1:nsub
            % megdata
            eval(['nam = ExpInfo.' vars{vvv} '.' datasets{ddd} '{sss};']);
            nam = [nam ext{vvv}];
            copyfile([rawpath in_files{vvv} filesep nam],[dspath out_files{vvv} filesep nam]);
        end
    end
end


%% copy the raw meg data to figshare
meg_path = [dspath 'MEGdata\'];
mkdir(dspath,'RawMEG')
cd(meg_path)
MyFolderInfo = dir(meg_path);
n_sub = length(MyFolderInfo)-2;

for i = 1:n_sub
   i
   filnam = MyFolderInfo(i+2).name;
   date_folder = dir(filnam);
   dat = date_folder(3).name;
   file_folder = dir([filnam '\' dat]);
   for f = 3:length(file_folder)
       fil = file_folder(f).name;
       if strcmp(fil(1),'b') %exclude data from the other two tasks
        movefile([filnam '\' dat '\' fil],[dspath 'RawMEG\20' dat '_' fil])
       end
   end
end

% delete the repetitive folder of epochs
cd(dspath)
rmdir(meg_path, 's')



%% get epoch data
epoch_path = [dspath 'epochs_tmp' filesep];
mkdir(epoch_path)
for sss = 1:nsub
    for ddd = 1:length(datasets)
        eval(['subnam = ExpInfo.subjects.' datasets{ddd} '{sss};']);
        mkdir(epoch_path, subnam)
        infile = ['Z:\Semantic\Analyse_data\' datasets{ddd} '_' subnam '\'];
        outfile = [epoch_path subnam filesep];
        copyfile([infile 'epoch_BL_Cross.mat'],[outfile 'epoch_BL_Cross.mat']);
        copyfile([infile 'epoch_WrdOn.mat'],[outfile 'epoch_WrdOn.mat']);
    end
end

%%copy the epoch data to figshare
addpath('Z:\fieldtrip-20220208\');
ft_defaults

mkdir(dspath,'epochs')
cd(epoch_path)
MyFolderInfo = dir(epoch_path);
n_sub = length(MyFolderInfo)-2;

for i = 1:n_sub
   i
   filnam = MyFolderInfo(i+2).name;
   load([epoch_path filnam '\epoch_WrdOn.mat']);
   valid_id = find(epoch.trialinfo(:,8) == 1); %only select the first fixation trials
   cfg = [];
   cfg.trials = valid_id;
   epoch = ft_selectdata(cfg, epoch);
   save([dspath 'epochs\' filnam '_WrdOn.mat'],'epoch')
   movefile([filnam '\epoch_BL_Cross.mat'],[dspath 'epochs\' filnam '_BL.mat'])  
end

% delete the repetitive folder of epochs
cd(dspath)
rmdir(epoch_path, 's')



% %% copy T1 images to figshare
% load('Z:\Lexical\Analyse_data\ExpInfo.mat')
% allsubs = ExpInfo.MEGcode_MRIcode(:,1);
% load('Z:\Lexical\data sharing\ExpInfo.mat')
% subs = ExpInfo.subjects.Targ60;
% valid_subid = cellfun(@(x) find(strcmp(allsubs, x)),subs);
% mri_id = ExpInfo1.MEGcode_MRIcode(valid_subid,:);
% % saveout
% ExpInfo.MEGcode_MRIcode = mri_id;
% save('Z:\Lexical\data sharing\ExpInfo.mat' ,'ExpInfo')
% % copy and rename mri images
% pdes= 'Z:\Lexical\data sharing\T1 images\';
% pmri = 'Z:\Lexical\RawData\MRI_images\';
% for i = 1:size(mri_id,1)
%     i
%     tmp = mri_id{i,2};
%     if ~strcmp(tmp,'nan')
%         fnii = 'T1_vol_v1_5.nii.gz';
%         pnii = [pmri tmp '_nifti\' tmp '_nifti\'];
%         if ~isfile([pnii fnii])
%             fnii(end-7) = '9';
%         end
%         copyfile([pnii fnii],pdes);
%         movefile([pdes fnii],[pdes tmp '_T1.nii.gz']);
%     end
% end


%% copy headmodels to figshare
mkdir(dspath,'HeadModel')
subs = ExpInfo.subjects.sv;
% copy headmodel to the sharing path
pdes = [dspath 'HeadModel\'];
pmri = 'Z:\Semantic\Results\Source_headmodel\';
for i = 1:size(subs,1)
    i
    copyfile([pmri 'Hdm_MRIaligned_' subs{i} '.mat'],pdes);
end





















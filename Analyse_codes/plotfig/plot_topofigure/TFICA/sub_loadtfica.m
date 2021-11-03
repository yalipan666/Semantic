function sub_loadtfica(hObject,eventdata,DirName,h_filename)
% load TF ICA and T-ICA data
global TF_IC T_IC PTF_IC PT_IC t f CommonChannelLabels pathname filename;

[filename,pathname] = uigetfile('*.mat','Select the MAT-file',DirName);
load ([pathname,filename])

TF_IC = [];
T_IC = [];
PTF_IC = [];
PT_IC = [];

filename = filename(1:end-4);
if exist('S_TF')
    TF_IC.P = (S_TF).*conj(S_TF)/Fs;%/U; %clear S_TF;
    TF_IC.s = s_TF;
    TF_IC.A = single(A_TF);
end
if exist('S_T')
    T_IC.P  = (S_T).*conj(S_T)/Fs;%/U; %clear S_T;
    T_IC.s = s_T;
    T_IC.A = single(A_T);
end
if exist('S_PTF')
    PTF_IC.P = (S_PTF).*conj(S_PTF)/Fs;%/U; %clear S_TF;
    PTF_IC.s = s_PTF;
    PTF_IC.A = single(A_PTF);
end
if exist('S_PT')
    PT_IC.P  = (S_PT).*conj(S_PT)/Fs;%/U; %clear S_T;
    PT_IC.s = s_PT;
    PT_IC.A = single(A_PT);
end

msgbox(['File ',filename,' is loaded sucessefully!'])
pause(1)
if ~isempty(h_filename)
set(h_filename,'String',filename);   
end
end
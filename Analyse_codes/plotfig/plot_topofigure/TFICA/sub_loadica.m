function sub_loadica(hObject,eventdata,DirName,h_filename)
% load ICA data
global T_IC PT_IC t f ChannelLabels pathname filename;

[filename,pathname] = uigetfile('*.mat','Select the MAT-file',DirName);
load ([pathname,filename])

T_IC = [];
PT_IC = [];

filename = filename(1:end-4);
if exist('S_T')
    T_IC.P  = (S_T).*conj(S_T)/Fs/U; %clear S_T;
    T_IC.s = s_T;
    T_IC.A = single(A_T);
end
if exist('S_PT')
    PT_IC.P  = (S_PT).*conj(S_PT)/Fs/U; %clear S_T;
    PT_IC.s = s_PT;
    PT_IC.A = single(A_PT);
end

msgbox(['File ',filename,' is loaded sucessefully!'])
pause(1)
if ~isempty(h_filename)
set(h_filename,'String',filename);   
end
end
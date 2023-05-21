%%
% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy,
%           Seoul National University, 2019.02.13.

%% Load files : Image files, Cell list, Blinder 1, Blinder 2.
clear
fnm='OE8';
load(['C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\' fnm '.mat'])
load(['C:\Users\Administrator\Dropbox\In vivo imaging\Behavior\' fnm '.mat'])

temp.Data=Final.Data;
temp.Deviation=Dev;
temp.Freeze=Freezing;
temp.Mouse=Final.Mouse;
temp.Mouse_info=Final.Mouse_info;
temp.behavior={'HC', 'CFC', 'Ret'};
clear Final
Final=temp;
save(['C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\New\' fnm '.mat'],'Final')
%%


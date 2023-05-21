% MODIFICATION HISTORY : 
%  2019.01.24.
%  Change the string finding condition while finding the mouse number in
%  excel file.
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.
%% Load data
clear
[fnm pth]= uigetfile('.mat','Pick the NF_map','multiselect','on');  % Cropped image folder
[numb txt raw]=xlsread('C:\Users\Administrator\Dropbox\In vivo imaging\Mouse_list'); %Mouse information
[fnm_fz pth_fz]= uigetfile('.mat','Behavior data','multiselect','on');  % Cropped image folder
Result_pth='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Classification\';
behavior={'HC','CFC','Ret','Ret2','Ret3'};%,'Ret4','Ret5','Ret6'};
overlapp=1; %1 for overlap 2 for not overlap.
%%
clear Data
if overlapp==1
for i=1:size(behavior,2)
    for j=1:length(fnm)
        sp_name=split(fnm{1,j},'_');
        spp_name=split(sp_name{end,1},'.');
    if ~isempty(strfind(char(spp_name{1,1}),char(behavior{1,i}))) && size(spp_name{1,1},2)==size(char(behavior{1,i}),2)
        load([pth char(fnm{1,j})])
    Data(:,1:3)=NF_map.list(:,1:3);
    Data(:,size(Data,2)+1)=NF_map.list(:,4);
    end
    end
end
else
    
   for i=1:size(behavior,2)
    for j=1:length(fnm)
    if ~isempty(strfind(char(fnm{1,j}),char(behavior{1,i}))) && ~isempty(strfind(char(behavior{1,i},char(fnm{1,j})))) % 2019.01.24
        load([pth char(fnm{1,j})])
    Data{i}(:,1:4)=NF_map.list(:,1:4);
    end
    end
   end 

end
%%
load([pth_fz fnm_fz])

%%
for j=1:size(raw,2)
        Mouse_information{1,j}=raw{1,j};
end
        
for i=1:size(raw,1)
    if ~isempty(strfind(Mouse_number,char(raw{i,1})))
        for j=1:size(raw,2)
        Mouse_information{2,j}=raw{i,j};
        end
    end
end
%%
Final.Data=Data;
Final.Distance=Distance;
Final.Freeze=Freezing;
Final.Mouse=Mouse_number;
Final.Mouse_info=Mouse_information;
Final.behavior=behavior;
save([Result_pth Mouse_number '.mat'],'Final')
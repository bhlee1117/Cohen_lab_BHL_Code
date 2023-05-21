%% Load tracking data
%%
%  Image enhancement code
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code traces mouse trajectory in squared shaped context.

% INPUTS 

% Image stack
% 

% OUTPUTS
% Filtered image stack

% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.


clear
[fnm pth]=uigetfile('.mat','Select the track data from CFC and Ret.','Multiselect','on');

bin=3;
CFC=[1 180 240];
Ret=[1 180];

%%
threshold=0.4;
    Freez_rate_Ret2=[];
    Freez_rate_Ret3=[];
    Freez_rate_Ret4=[];
    Freez_rate_Ret5=[];
    Freez_rate_Ret6=[];
for i=1:length(fnm)
    
    data=importdata([pth char(fnm{1,i})]);
    track{i}=data.track;
    tmp=split(char(fnm{1,i}),'.');
    spname=split(tmp{1,1},'_');
    behave=char(spname{2,1});

    switch behave
        case 'Ret'
    dist{i}=data.ave_disp;        
    clear desk
    desk=data.ave_disp(ceil(Ret(1,1)/bin):Ret(1,2)/bin,:);
    Freez_rate_Ret=sum(desk<threshold)/size(desk,1);
    
        case 'CFC'
            dist{i}=data.ave_disp;        
    clear desk
    desk=data.ave_disp(ceil(CFC(1,1)/bin):CFC(1,2)/bin,:);
    Freez_rate_BS=sum(desk<threshold)/size(desk,1);
    clear desk
    desk=data.ave_disp(ceil(CFC(1,2)/bin):CFC(1,3)/bin,:);
    Freez_rate_CFC=sum(desk<threshold)/size(desk,1);
    
     case 'Ret2'
            dist{i}=data.ave_disp;        
    dist{i}=data.ave_disp;        
    clear desk
    desk=data.ave_disp(ceil(Ret(1,1)/bin):Ret(1,2)/bin,:);
    Freez_rate_Ret2=sum(desk<threshold)/size(desk,1);
   
    
    case 'Ret3'
            dist{i}=data.ave_disp;        
    dist{i}=data.ave_disp;        
    clear desk
    desk=data.ave_disp(ceil(Ret(1,1)/bin):Ret(1,2)/bin,:);
    Freez_rate_Ret3=sum(desk<threshold)/size(desk,1);
    
    case 'Rem'
            dist{i}=data.ave_disp;        
    dist{i}=data.ave_disp;        
    clear desk
    desk=data.ave_disp(ceil(Ret(1,1)/bin):Ret(1,2)/bin,:);
    Freez_rate_Ret4=sum(desk<threshold)/size(desk,1);
    
    case 'Ret5'
            dist{i}=data.ave_disp;        
    dist{i}=data.ave_disp;        
    clear desk
    desk=data.ave_disp(ceil(Ret(1,1)/bin):Ret(1,2)/bin,:);
    Freez_rate_Ret5=sum(desk<threshold)/size(desk,1);
    
    case 'Ret6'
            dist{i}=data.ave_disp;        
    dist{i}=data.ave_disp;        
    clear desk
    desk=data.ave_disp(ceil(Ret(1,1)/bin):Ret(1,2)/bin,:);
    Freez_rate_Ret6=sum(desk<threshold)/size(desk,1);
    end
end
%%
Freezing.Baseline=Freez_rate_BS;
Freezing.CFC=Freez_rate_CFC;
Freezing.Retrieval=Freez_rate_Ret;
if ~isempty(Freez_rate_Ret2)
    Freezing.Retrieval2=Freez_rate_Ret2;
    Distance.Retrieval2=dist{1,3};
end
if ~isempty(Freez_rate_Ret3)
    Freezing.Retrieval3=Freez_rate_Ret3;
    Distance.Retrieval3=dist{1,4};
end
if ~isempty(Freez_rate_Ret4)
    Freezing.Retrieval4=Freez_rate_Ret4;
    Distance.Retrieval4=dist{1,5};
end

if ~isempty(Freez_rate_Ret5)
    Freezing.Retrieval5=Freez_rate_Ret5;
    Distance.Retrieval5=dist{1,6};
end
if ~isempty(Freez_rate_Ret6)
    Freezing.Retrieval6=Freez_rate_Ret6;
    Distance.Retrieval6=dist{1,7};
end
Distance.CFC=dist{1,1};
Distance.Retrieval=dist{1,2};
Distance.track=track;
Mouse_number=spname{1,1};
save([pth char(spname{1,1}) '.mat'],'Freezing','Distance','Mouse_number')
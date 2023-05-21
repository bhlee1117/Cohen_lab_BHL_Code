%% Load tracking data
%%
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code evaluates the freezing rate of the mouse after behavior
% tracking.

% INPUTS 

% Dev_track '*.mat' files
% 

% OUTPUTS
% Filtered image stack

% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.


clear
[fnm pth]=uigetfile('.mat','Select the track data from CFC and Ret.','Multiselect','on');

CFC=[1 180 240];
Ret=[1 180];
fps=4;
bout_time=0.5;%sec
%%
threshold=3.7;
Freez_rate_BS=[];
Freez_rate_CFC=[];
    Freez_rate_Ret2=[];
    Freez_rate_Ret3=[];
    Freez_rate_Ret4=[];
    Freez_rate_Ret5=[];
    Freez_rate_Ret6=[];
    Freez_rate_Rem=[];
    bout_n=fps*bout_time;
for i=1:length(fnm)
    
    data=importdata([pth char(fnm{1,i})]);
    tmp=split(char(fnm{1,i}),'.');
    spname=split(tmp{1,1},'_');
    behave=char(spname{2,1});

    switch behave
        case 'Ret1'
            Dev{i}=data;
            clear fz_length
    fz_length=bin_length(data<threshold);
    Freez_rate_Ret=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/(Ret(1,2)*fps);
    
        case 'CFC'
            clear fz_length
               Dev{i}=data;
    fz_length=bin_length(data(ceil(CFC(1,1)*fps):CFC(1,2)*fps,1)<threshold);
    Freez_rate_BS=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/((CFC(1,2)-CFC(1,1)+1)*fps);
    clear fz_length
    fz_length=bin_length(data(ceil(CFC(1,2)*fps):CFC(1,3)*fps,1)<threshold);
    Freez_rate_CFC=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/((CFC(1,3)-CFC(1,2)+1)*fps);

        case 'Ret2'
            Dev{i}=data;
            clear fz_length
    fz_length=bin_length(data<threshold);
    Freez_rate_Ret2=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/(Ret(1,2)*fps);
   
    
    case 'Ret3'
            Dev{i}=data;
            clear fz_length
    fz_length=bin_length(data<threshold);
    Freez_rate_Ret3=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/(Ret(1,2)*fps);
    
    case 'Ret4'
            Dev{i}=data;
            clear fz_length
    fz_length=bin_length(data<threshold);
    Freez_rate_Ret4=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/(Ret(1,2)*fps);
    
    case 'Rem'
           Dev{i}=data;
            clear fz_length
    fz_length=bin_length(data<threshold);
    Freez_rate_Rem=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/(Ret(1,2)*fps);
    
    case 'Ret5'
            Dev{i}=data;
            clear fz_length
    fz_length=bin_length(data<threshold);
    Freez_rate_Ret5=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/(Ret(1,2)*fps);
    
    case 'Ret6'
            Dev{i}=data;
            clear fz_length
    fz_length=bin_length(data<threshold);
    Freez_rate_Ret6=sum(fz_length(find(fz_length(:,2)>=bout_n),2))/(Ret(1,2)*fps);
    
    end
end



Freezing.Baseline=Freez_rate_BS;
Freezing.CFC=Freez_rate_CFC;
Freezing.Retrieval=Freez_rate_Ret;
if ~isempty(Freez_rate_Ret2)
    Freezing.Retrieval2=Freez_rate_Ret2;
end
if ~isempty(Freez_rate_Ret3)
    Freezing.Retrieval3=Freez_rate_Ret3;

end
if ~isempty(Freez_rate_Ret4)
    Freezing.Retrieval4=Freez_rate_Ret4;
    
end

if ~isempty(Freez_rate_Ret5)
    Freezing.Retrieval5=Freez_rate_Ret5;
    
end
if ~isempty(Freez_rate_Ret6)
    Freezing.Retrieval6=Freez_rate_Ret6;
  
end

if ~isempty(Freez_rate_Rem)
    Freezing.Remote=Freez_rate_Rem;
  
end

save([pth char(spname{1,1}) '.mat'],'Freezing','Dev')


% Distance.CFC=dist{1,1};
% Distance.Retrieval=dist{1,2};
% Distance.track=track;
% Mouse_number=spname{1,1};
% save([pth char(spname{1,1}) '.mat'],'Freezing','Distance','Mouse_number')
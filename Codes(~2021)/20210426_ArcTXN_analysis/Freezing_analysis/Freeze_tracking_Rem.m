%% Load tracking data
clear
[fnm pth]=uigetfile('.mat','Select the track data from CFC, Ret.& Remote','Multiselect','on');

bin=3;
CFC=[1 180 240];
Ret=[1 180];
Rem=[1 180];
%%
threshold=18;
for i=1:length(fnm)
    
    data=importdata([pth char(fnm{1,i})]);
    
    tmp=split(char(fnm{1,i}),'.');
    spname=split(tmp{1,1},'_');
    behave=char(spname{2,1});
    
    switch behave
        case 'Rem'
    dist{3}=data.ave_disp;     
    track{3}=data.track;
    clear desk
    desk=data.ave_disp(ceil(Rem(1,1)/bin):Rem(1,2)/bin,:);
    Freez_rate_Rem=sum(desk<threshold)/size(desk,1);
    
    case 'Ret'
    dist{2}=data.ave_disp;
    track{2}=data.track;
    clear desk
    desk=data.ave_disp(ceil(Ret(1,1)/bin):Ret(1,2)/bin,:);
    Freez_rate_Ret=sum(desk<threshold)/size(desk,1);
    
        case 'CFC'
            dist{1}=data.ave_disp;   
            track{1}=data.track;
    clear desk
    desk=data.ave_disp(ceil(CFC(1,1)/bin):CFC(1,2)/bin,:);
    Freez_rate_BS=sum(desk<threshold)/size(desk,1);
    clear desk
    desk=data.ave_disp(ceil(CFC(1,2)/bin):CFC(1,3)/bin,:);
    Freez_rate_CFC=sum(desk<threshold)/size(desk,1);
    
    end
end
%%
Freezing.Baseline=Freez_rate_BS;
Freezing.CFC=Freez_rate_CFC;
Freezing.Retrieval=Freez_rate_Ret;
Freezing.Remote=Freez_rate_Rem;
Distance.CFC=dist{1,1};
Distance.Retrieval=dist{1,2};
Distance.Remote=dist{1,3};
Distance.track=track;
Mouse_number=spname{1,1};
save([pth char(spname{1,1}) '_w_rem.mat'],'Freezing','Distance','Mouse_number')
function [Conv_trac A_signal_trace]=activitiy_generator(N,basal_level)
%% set parameters
time_bin=1; %min
period=60*50; %min
t=[time_bin:time_bin:period];
iter=10;
CFC_response_rate=25;

peak_GFP=100; decay_GFP=120; shift_GFP=20; end_time_GFP=14*60;
peak_TXN=6; decay_TXN=2.5; shift_TXN=0; end_time_TXN=20;
CFC_time= [60*24-3:60*24];

%%
clear Conv_trac   A_signal_trace

A_signal_trace=double(rand(N,length(t))>(1-basal_level/100));
if ~isempty(CFC_time)
    CFC_rand=rand(N,1);
    [p order]=sort(CFC_rand,'ascend');
    CFC_pop=CFC_rand<CFC_rand(order(round(N/100*CFC_response_rate)));
    %CFC_signal=reshape(CFC_rand<CFC_rand(order(round(N/100*CFC_response_rate))),[],4);
    %size(order(1:round(N/100*CFC_response_rate))')
    CFC_signal=zeros(N,size(CFC_time,2));
    for j=find(CFC_pop)'
        CFC_signal(j,ceil(rand*size(CFC_time,2)))=1;
    end
    A_signal_trace(:,CFC_time/time_bin)=CFC_signal;
end

for res=1:2
Conv_trac{1,res}=zeros(size(A_signal_trace,1),size(A_signal_trace,2));
end
%reponse_ftn{res}=reponse_ftn{res}/max(reponse_ftn{res});

for n=1:size(A_signal_trace,1)
    S=find(A_signal_trace(n,:));
    
    for ss=S
    tmp=zeros(1,size(A_signal_trace,2));    
        tmp(ss)=1;
        reponse_ftn{1}=generate_response_ftn(peak_GFP+randn*peak_GFP*0.1,decay_GFP+randn*30,end_time_GFP,shift_GFP);
        reponse_ftn{2}=generate_response_ftn(peak_TXN+randn,decay_TXN+randn,end_time_TXN,shift_TXN);
        for res=1:2
            reponse_ftn{res}=[zeros(1,size(reponse_ftn{res},2)) reponse_ftn{res}];
            reponse_ftn{res}(reponse_ftn{res}<0)=0;
            Conv_trac{1,res}(n,:)=Conv_trac{1,res}(n,:)+conv2(tmp,reponse_ftn{res},'same');    
        end
    end
    
end

end
%end
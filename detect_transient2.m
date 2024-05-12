function [trans transients_final]=detect_transient2(traces_hi,thres,spike_trace,Prox_spike_tau)

thres_on=thres(1);
thres_off=thres(2);
if nargin<4
Prox_spike_tau=15; %how close spike will be combined to the transient
end
for i=1:size(traces_hi,1) %neuron

    [transients(i,:) n]=bwlabel(traces_hi(i,:)>thres_on);
    [transients_back(i,:)]=bwlabel(traces_hi(i,:)>thres_off);
    transients_final(i,:)=double(movmean(traces_hi(i,:),3)>thres_on);
    %transients_final(i,:)=zeros(1,size(traces_hi,2));
    tr_list=[1:n]; rm=[];
    sp_time=find(spike_trace(i,:));
    for t=1:n %transients
        tr_ind=find(transients(i,:)==t);
        t2=find(transients_back(i,:)==transients_back(i,tr_ind(end)));
        %transients_final(tr_ind(1):t2(end))=1;
        sp_time_trans=find(spike_trace(i,t2(1):t2(end)));
        if ~isempty(sp_time_trans)
            nearest_ISI=(sp_time_trans(1)+t2(1)-1)-sp_time;
            include_spike=find(nearest_ISI>0 & nearest_ISI<Prox_spike_tau);
        while ~isempty(include_spike)
            t2(1)=sp_time(include_spike(1));
            sp_time_trans=find(spike_trace(i,t2(1):t2(end)));
            nearest_ISI=(sp_time_trans(1)+t2(1)-1)-sp_time;
            include_spike=find(nearest_ISI>0 & nearest_ISI<Prox_spike_tau);
        end
        end
        transients_final(i,t2(1):t2(end))=1;
    end
    [transients_final(i,:) n]=bwlabel(double(transients_final(i,:)));
    if n>1
        for t=1:n %transients
            tr_ind=find(transients_final(i,:)==t);
            trans(i).length(t)=length(tr_ind);
            trans(i).amp(t)=max(traces_hi(i,tr_ind));
            trans(i).int(t)=sum(traces_hi(i,tr_ind));
            sp_seg=spike_trace(i,tr_ind);
            trans(i).spike_number(t)=sum(sp_seg);
            if sum(sp_seg)>1
                sp_seg_time=find(sp_seg);
                trans(i).mean_ISI(t)=mean(sp_seg_time(2:end)-sp_seg_time(1:end-1));
                trans(i).ISI{t}=sp_seg_time(2:end)-sp_seg_time(1:end-1);
            else
                trans(i).mean_ISI(t)=NaN;
                trans(i).ISI{t}=NaN;
            end
        end
        trans(i).interval=transients_final(i,:);
    else
        trans(i).length=[];
        trans(i).amp=[];
        trans(i).int=[];
        trans(i).spike_number=[];
        trans(i).mean_ISI=[];
        trans(i).interval=[];
    end
end
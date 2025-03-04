function [trans transients_final]=detect_transient(traces_hi,thres,spike_trace)

thres_on=thres(1);
thres_off=thres(2);

for i=1:size(traces_hi,1) %neuron
    [transients(i,:) n]=bwlabel(traces_hi(i,:)>thres_on);
    [transients_back(i,:)]=bwlabel(traces_hi(i,:)>thres_off);
    transients_final(i,:)=double(movmean(traces_hi(i,:),3)>thres_on);
    tr_list=[1:n]; rm=[];
    for t=1:n %transients
        tr_ind=find(transients(i,:)==t);
        t2=find(transients_back(i,:)==transients_back(i,tr_ind(end)));
        %transients_final(tr_ind(1):t2(end))=1;
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
            else
                trans(i).mean_ISI(t)=NaN;
            end
        end
        trans(i).interval=transients_final(i,:);
    else
        trans(i).length=[];
        trans(i).amp=[];
        trans(i).int=[];
        trans(i).spike_number=[];
        trans(i).mean_ISI=[];
        trans(i).interval=zeros(1,size(traces_hi,2));
    end
end
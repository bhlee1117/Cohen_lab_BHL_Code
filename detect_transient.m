function [trans tr_trace]=detect_transient(traces_hi,burst_pow)
tr_trace=[];
for i=1:size(traces_hi,1) %neuron
    i
    thres=get_threshold(traces_hi(i,:),0.9); thres2=thres*1/2;
    [transients n]=bwlabel(movmean(traces_hi(i,:),7)>thres);
    [transients_back]=bwlabel(movmean(traces_hi(i,:),7)>thres2);
    transients_final=movmean(traces_hi(i,:),3)>thres;
    tr_list=[1:n]; rm=[];
    for t=1:n
    tr_ind=find(transients==t);
    t2=find(transients_back==transients_back(tr_ind(end)));
    transients_final(tr_ind(1):t2(end))=1;
    end
    [transients_final n]=bwlabel(transients_final);
    
    for t=1:n
    tr_ind=find(transients_final==t);
    
    trans(i).length(t)=length(tr_ind);
    trans(i).amp(t)=max(traces_hi(i,tr_ind));
    trans(i).int(t)=sum(traces_hi(i,tr_ind));
    trans(i).burst(t)=sum(burst_pow(i,tr_ind));
    end
    trans(i).interval=transients_final;
    tr_trace=[tr_trace; transients_final];
end
end
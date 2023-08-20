%% Complex spike -> single spike adaptation figure generation

% load data
clear
cd('/Users/bhlee1117/Dropbox/20230718')
Voltron1_result=importdata("20230521_Adaptation_result_BHLm2526.mat");
Voltron2_result=importdata("Result_V2CheRiffST_20230714.mat");

%% Scratch good cells

V2=Voltron2_result.Result;
V2(Voltron2_result.l)=[];
V2Path=Voltron2_result.fpath;
V2Path(Voltron2_result.l)=[];

for i=1:size(V2,2)
    sp=split(V2Path{i},'_');
    cstr=find(contains(sp,'ell'));
    pathstr{1,i}=cell2mat(sp(cstr-1:cstr)');
end
[pathunq patharg] = unique(pathstr, 'stable');

for i=1:length(pathunq)

    session = find(strcmp(pathstr, pathunq{i}));
    BlueStimNumber=[];
    for s=1:length(session)
        V2Result(i,s)=V2(session(s));
        V2Result{i,s}.fpath=V2Path{session(s)};
        BlueStimNumber(s)=max(bwlabel(V2{session(s)}.Blue>0));
    end
    [~, arg]=sort(BlueStimNumber,'descend');
    V2Result(i,1:length(session))=V2Result(i,arg);
end

%%
figure(1); clf;
scale=20;
ax1=nexttile([1 1]);
clear Transients_interval Transients CS_wvform CS_spike
CS_wvform=[]; CS_blue=[]; CS_spike=[];
for i=1:size(V2Result,1)
Sessions=find(1-cellfun(@isempty,V2Result(i,:)));
    tr=V2Result{i,Sessions(end)}.trace;
    tr=tr-median(tr); tr=tr./get_threshold(tr,1);
    tr_CS=tr-movmedian(tr,300,2); % high pass
    tr_CS=tr_CS/get_threshold(tr_CS,1); % normalize to noise level

    tr_SS=tr-movmedian(tr,30,2); % high pass
    tr_SS=tr_SS/get_threshold(tr_SS,1); % normalize to noise level

            spike_SS=find_spike_bh(tr_SS,4,3);
        spike_complex=find_spike_bh(tr_CS,3,1);

    [trans tr_trace]=detect_transient(movmean(tr_CS,5,2),[3 1]);
    [CS_list CS]=find_CS(tr_CS,spike_SS,11,2); % t_threshold area
    ComplexTransient=bwlabel((CS+tr_trace)>0);
    Transients=zeros(max(ComplexTransient),4);
    for t=1:max(ComplexTransient)
        tr_ind=find(ComplexTransient==t);
        Transients(t,:)=[length(tr_ind) max(tr_CS(1,tr_ind)) sum(tr_CS(1,tr_ind)) sum(spike_SS(tr_ind))];
    end

    rmv_trans=find(Transients(:,1)<20 | Transients(:,3)<15 | Transients(:,4)<2);
    ComplexTransient(ismember(ComplexTransient,rmv_trans))=0;
    ComplexTransient=bwlabel(ComplexTransient>0);

    spike_complex(ComplexTransient==0)=0;
    spike_total=spike_SS | spike_complex;

    for t=1:max(ComplexTransient)
        wvform=zeros(1,120);
        sp_tmp=NaN(1,20);
        tr_ind=find(ComplexTransient==t);
        sp_ind=find(spike_total(tr_ind));
        Transients_interval(t,:)=[tr_ind(1) tr_ind(end) tr_ind(sp_ind(1)) sum(sp_ind)]; %transiend start, end, first spike
        tmp=tr(Transients_interval(t,3)-10:tr_ind(end));
        wvform(1:length(tmp))=tmp;
        CS_wvform=[CS_wvform; wvform];
        CS_blue=[CS_blue; V2Result{i,Sessions(end)}.Blue(Transients_interval(t,3))];
        sp_tmp(1:sum(spike_total(tr_ind)))=[sum(spike_total(tr_ind)) sp_ind(2:end)-sp_ind(1:end-1)];
        CS_spike=[CS_spike; sp_tmp];
    end


        tr_CS=tr; tr_CS(~(spike_total & ComplexTransient>0))=NaN;
        tr_SS=tr; tr_SS(~(spike_total & ComplexTransient==0))=NaN;

        plot_label([1:length(tr)],tr+scale*i,ComplexTransient>0,[0 0 0])
        hold all

        plot([1:length(tr)],tr_CS+scale*i,'marker','.','color',[0.1 0.7 1])
        plot([1:length(tr)],tr_SS+scale*i,'marker','.','color',[1 0 0])

% Wavelet transform
Fs = 1000; % Sampling frequency (change as per your signal)
t = 0:1/Fs:1; % Time vector

waveletName = 'cmor'; % Wavelet family (you can choose other wavelets as well)
scales = 1:128; % Range of scales (adjust as needed)

[coefficients{i}, frq{i}] = cwt(tr, 1000);

end

%%


%%
i=5;
ax2=nexttile([1 1]);
imagesc(abs(coefficients{i}));

set(gca, 'ytick', [1:5:size(coefficients{i},1)], 'yticklabel' , round(frq{i}(1:5:end))');
colormap('jet');
colorbar;
xlabel('Time (s)');
ylabel('Scale');
linkaxes([ax1 ax2],'x')


%% Cs types
[blue_sorted blue_arg]=sort(CS_blue,'ascend');
figure;
s=[1:70:548];
cmap=jet(length(s));
for i=1:length(s)-1
plot(mean(CS_wvform(blue_arg(s(i):s(i+1)),:),1),'color',cmap(i,:))
hold all
end

figure;
nexttile([1 1])
plot(blue_sorted,CS_spike(blue_arg,1),'.')

nexttile([1 1])
cmap=jet(5);
for i=1:5

plot(blue_sorted,CS_spike(blue_arg,i+1),'o','color',cmap(i,:))    
hold all
end

nexttile([1 1])
cmap=distinguishable_colors(5);
for i=1:5
    [h b]=histcounts(CS_spike(~isnan(CS_spike(:,i+1)),i+1),[1:1:15],'Normalization','probability');
    plot(b(2:end)-0.5,h,'color',cmap(i,:),'marker','.')
    hold all
end

nexttile([1 1])
cmap=distinguishable_colors(5);
blue_bin=[0:max(blue_sorted)/5:max(blue_sorted)];
for i=1:length(blue_bin)-1
    find(blue_bin(i)>CS_blue & blue_bin(i)<CS_blue)
    CS_spike(find(blue_bin(i)>CS_blue & blue_bin(i)<CS_blue),2)
    plot(b(2:end)-0.5,h,'color',cmap(i,:),'marker','.')
    hold all
end

 figure;
for i=1:4
nexttile([1 1])
plot(CS_blue,CS_spike(:,i+1),'.','color',cmap(i,:))
tmp=~isnan(CS_spike(:,i+1));
[r p]=corr(CS_blue(tmp),CS_spike(tmp,i+1));
title(sprintf('Correlation = %.2f, p = %.4f', r, p))
hold all
xlabel('Blue intensity (V)')
ylabel('ISI (ms)')
end
%%
for i=1:size(NiceNeuronList,2)
    ax1=[];
    %ax1=[ax1 nexttile([1 1])];
    Av=find(AvailResult(NiceNeuronList(1,i),:));

    for j=find(Waveform_type(NiceNeuronList(1,i),:)==1)%[NiceNeuronList(2,i) setdiff(Av,NiceNeuronList(2,i))]
        t=[1:length(Result{NiceNeuronList(1,i),j}.blue)];
        tr=Result{NiceNeuronList(1,i),j}.traces(t);
        tr=tr-median(tr);

        tr=tr/get_threshold(tr,1);
        sp=find(Result{NiceNeuronList(1,i),j}.spike(t));
        tr_transient_hi=tr-movmedian(tr,300,2);
        [trans tr_trace]=detect_transient(movmean(tr_transient_hi,20,2),[2 1]);
        [CS_list CS]=find_CS(tr-movmedian(tr,150,2),Result{NiceNeuronList(1,i),j}.spike(t),15,5);
        Result{NiceNeuronList(1,i),j}.Complex=CS>0;
        CS=CS>0; tr_CS=tr; tr_CS(~CS)=NaN;
        tr_Trns=tr; tr_Trns(tr_trace==0)=NaN;
        plot(t,tr+g*20,'color',cmap(i,:))
        hold all
        plot(sp,tr(sp)+g*20,'r.')
        plot(t,tr_Trns+g*20,'color',[1 0 0.7])
        plot(t,tr_CS+g*20,'color',[1 0 0.7])


        g=g+1;
        % ax1=[ax1 nexttile([1 1])];
        % plot(Result{NiceNeuronList(1,i),j}.blue)
        % linkaxes([ax1],'x')
    end
end
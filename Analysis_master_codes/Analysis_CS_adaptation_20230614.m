
clear
cd /Users/bhlee1117/Documents/BHL/Matlab_project/20230614_Adaptation
load('20230521_Adaptation_result_BHLm2526.mat')
% NiceNeuronList=sum(1-cellfun(@isempty,Result),2);
% list_pyr=find(NiceNeuronList>1)';
%y=[8 16 19 23 27 29 33 34 36 37 41 42 45 46 48 60 61 65 68 69 70 73 83 84 88 90 97 102 104 107 110 114 117 118];
Blue_CAL=[[0.01:0.005:5]' interp1(Blue_CAL(:,1),Blue_CAL(:,2),[0.01:0.005:5])'];
%% Show cell images
figure;
for i=NiceNeuronList(1,:)
    nexttile([1 1])
    imshow2(Result{i,1}.ref_im,[])
    nexttile([1 1])
    imshow2(Result{i,1}.c_ftprnt(:,:,1),[])
end
Waveform_type=zeros(size(Result)); % 1=Ramp pulse, 2= 3 pulses
for i=1:size(Result,1)
    for j=1:size(Result,2)
    if ~isempty(Result{i,j})
        if max(bwlabel(Result{i,j}.blue))>3
        Waveform_type(i,j)=1;    
        else
        Waveform_type(i,j)=2;    
        end
    else
    Waveform_type(i,j)=NaN;
    end
    end
end
%%
AvailResult=1-cellfun(@isempty,Result);
cmap=winter(size(NiceNeuronList,2));
g=1;
figure;
tick=[];

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


%% Ramp region

figure;
g=1; scale=17;

Transient_traces=[]; FineSpikes=[];
for i=1:23%size(NiceNeuronList,2)
%for i=24:size(NiceNeuronList,2)
    ax1=[];
    %ax1=[ax1 nexttile([1 1])];

    for j=find(Waveform_type(NiceNeuronList(1,i),:)==1)
        Blue_sessions=bwlabel(Result{NiceNeuronList(1,i),j}.blue);
        t=find(Blue_sessions==max(Blue_sessions)); % Ramp time
        tr=Result{NiceNeuronList(1,i),j}.traces(t);
        tr=tr-median(tr); tr=tr./get_threshold(tr,1);
        
        tr_hi=tr-movmedian(tr,300,2); % high pass
        tr_hi=tr_hi/get_threshold(tr_hi,1); % normalize to noise level
        tr_hi2=tr-movmedian(tr,30,2); % high pass
        tr_hi2=tr_hi2/get_threshold(tr_hi2,1); % normalize to noise level
        spike_in=find_spike_bh(tr_hi2,4,3);
        spike_complex=find_spike_bh(tr_hi,3,2);
        [trans tr_trace]=detect_transient(movmean(tr_hi,5,2),[3 1]);
        [CS_list CS]=find_CS(tr_hi,spike_in,15,1); % t_threshold area
        ComplexTransient=bwlabel((CS+tr_trace)>0);

        Transients=zeros(max(ComplexTransient),4);
        for t=1:max(ComplexTransient)
            tr_ind=find(ComplexTransient==t);
            Transients(t,:)=[length(tr_ind) max(tr_hi(1,tr_ind)) sum(tr_hi(1,tr_ind)) sum(spike_in(tr_ind))];
        end

        rmv_trans=find(Transients(:,1)<20 | Transients(:,3)<15 | Transients(:,4)<2);
        ComplexTransient(ismember(ComplexTransient,rmv_trans))=0;
        ComplexTransient=bwlabel(ComplexTransient>0);
        spike_complex(ComplexTransient==0)=0;
        spike_total=spike_in | spike_complex;

        tr_CS=tr; tr_CS(~(spike_total & ComplexTransient>0))=NaN;
        tr_SS=tr; tr_SS(~(spike_total & ComplexTransient==0))=NaN;

        plot_label([1:length(tr)],tr+scale*g,ComplexTransient>0,[0 0 0])
        hold all

        plot([1:length(tr)],tr_CS+scale*g,'marker','.','color',[0.1 0.7 1])
        plot([1:length(tr)],tr_SS+scale*g,'marker','.','color',[1 0 0])
        
        Transient_traces(g,:)=ComplexTransient;
        FineSpikes(g,:)=spike_total;
        %Result{i,j}.FineSpikes=spike_total;

        g=g+1;
    end     
end

figure; 
bin_size=150;
clf; cmap=distinguishable_colors(2);
CS_trace=FineSpikes; CS_trace(Transient_traces==0)=0;
SS_trace=FineSpikes; SS_trace(Transient_traces~=0)=0;
CS=movsum(CS_trace,bin_size,2); CS=CS(:,[1:bin_size:end]);%./sum(FineSpikes,2);
SS=movsum(SS_trace,bin_size,2);SS=SS(:,[1:bin_size:end]);%./sum(FineSpikes,2);
N=size(FineSpikes,1);
errorbar([1:bin_size:size(Transient_traces,2)]/1000,mean(CS,1),std(CS,0,1)/sqrt(N),'marker','+','color',cmap(1,:))
hold all
errorbar([1:bin_size:size(Transient_traces,2)]/1000,mean(SS,1),std(SS,0,1)/sqrt(N),'marker','+','color',cmap(2,:))
legend('Complex spike','Single spike')
xlabel('Time (sec)')
%ylabel('Fraction of spike')
ylabel('Number of spike')

%% Pulse region
figure;
g=1; scale=17;
Transient_traces=[]; FineSpikes=[]; NormBlue=[];

for i=1:size(NiceNeuronList,2)
%for i=24:size(NiceNeuronList,2)
    ax1=[];
    %ax1=[ax1 nexttile([1 1])];

    for j=find(Waveform_type(NiceNeuronList(1,i),:)==1) % sessions
        Blue_sessions=bwlabel(Result{NiceNeuronList(1,i),j}.blue); % Blue stim regions
        
        % Find first spiking blue intensity
        delay=20; gg=1; b=1; % Spike showing in 20 ms after stim
        while gg
            blue_seg=find(Blue_sessions==b);
            if sum(Result{NiceNeuronList(1,i),j}.spike(blue_seg(1):blue_seg(end)+delay))>0
                gg=0;
                BlueVoltage=Result{NiceNeuronList(1,i),j}.blue(blue_seg(1));
                FirstSpikeBlue(g,1)=Blue_CAL(round(Blue_CAL(:,1)*1000)==round(BlueVoltage*1000),2);
            end
            b=b+1;
        end

        for Pulse=[5:-1:1]
        BlueTime=find(Blue_sessions==max(Blue_sessions)-Pulse); % Ramp time
        BlueTime=[BlueTime(1):BlueTime(1)+500]; % First 500 ms
        tr=Result{NiceNeuronList(1,i),j}.traces(BlueTime);
        tr=tr-median(tr); tr=tr./get_threshold(tr,1);
        
        tr_hi=tr-movmedian(tr,300,2); % high pass
        tr_hi=tr_hi/get_threshold(tr_hi,1); % normalize to noise level
        tr_hi2=tr-movmedian(tr,30,2); % high pass
        tr_hi2=tr_hi2/get_threshold(tr_hi2,1); % normalize to noise level
        spike_in=find_spike_bh(tr_hi2,4,3);
        spike_complex=find_spike_bh(tr_hi,3,2); % Find with lower threshold (for complex spikes)
        [trans tr_trace]=detect_transient(movmean(tr_hi,5,2),[3 1]);
        [CS_list CS]=find_CS(tr_hi,spike_in,15,1); % t_threshold area
        ComplexTransient=bwlabel((CS+tr_trace)>0);

        Transients=zeros(max(ComplexTransient),4);
        for t=1:max(ComplexTransient)
            tr_ind=find(ComplexTransient==t);
            Transients(t,:)=[length(tr_ind) max(tr_hi(1,tr_ind)) sum(tr_hi(1,tr_ind)) sum(spike_in(tr_ind))];
        end

        rmv_trans=find(Transients(:,1)<20 | Transients(:,3)<15 | Transients(:,4)<2);
        ComplexTransient(ismember(ComplexTransient,rmv_trans))=0;
        ComplexTransient=bwlabel(ComplexTransient>0);
        spike_complex(ComplexTransient==0)=0;
        spike_total=spike_in | spike_complex;

        tr_CS=tr; tr_CS(~(spike_total & ComplexTransient>0))=NaN;
        tr_SS=tr; tr_SS(~(spike_total & ComplexTransient==0))=NaN;

        subplot(1,5,6-Pulse)
        plot_label([1:length(tr)],tr+scale*g,ComplexTransient>0,[0 0 0])
        hold all
        plot([1:length(tr)],tr_CS+scale*g,'marker','.','color',[0.1 0.7 1])
        plot([1:length(tr)],tr_SS+scale*g,'marker','.','color',[1 0 0])

        Transient_traces(g,:,6-Pulse)=ComplexTransient;
        FineSpikes(g,:,6-Pulse)=spike_total;
        BlueVoltage=Result{NiceNeuronList(1,i),j}.blue(BlueTime(1)); BlueIntens=Blue_CAL(find(round(Blue_CAL(:,1)*1000)==round(BlueVoltage*1000)),2);
        NormBlue(g,6-Pulse)=BlueIntens/FirstSpikeBlue(g,1);
        end

        g=g+1;
    end     
end

Blue_bin=[0 3 6 10 15 30]; cmap=jet(length(Blue_bin)-1);
figure; 
bin_size=100;
clear CS_trace SS_trace totalSpike
for b=2:length(Blue_bin)

[targetN targetP]=find(NormBlue>Blue_bin(b-1) & NormBlue<=Blue_bin(b));
if ~isempty(targetN)
for n=1:size(targetN,1)
CS_trace{b}(n,:)=FineSpikes(targetN(n),:,targetP(n));
CS_trace{b}(n,Transient_traces(targetN(n),:,targetP(n))==0)=0;

SS_trace{b}(n,:)=FineSpikes(targetN(n),:,targetP(n));
SS_trace{b}(n,Transient_traces(targetN(n),:,targetP(n))~=0)=0;

totalSpike{b}(n,:)=sum(FineSpikes(targetN(n),:,targetP(n)));
end
ax1=subplot(1,4,1);
title('Complex spikes')
CS=movsum(CS_trace{b},bin_size,2); CS=CS(:,[1:bin_size:end]);
%plot([1:bin_size:size(FineSpikes,2)]/1000,CS,'color',cmap(b-1,:))
%errorbar([1:bin_size:size(FineSpikes,2)]/1000,mean(CS,1),std(CS,0,1)/sqrt(size(CS,1)),'marker','+','color',cmap(b-1,:))
plot([1:bin_size:size(FineSpikes,2)]/1000,mean(CS,1,'omitnan'),'marker','*','color',cmap(b-1,:))
hold all
xlabel('Time (sec)')
ylabel('Number of spike')

ax3=subplot(1,4,3);
title('Complex spikes')
CS=movsum(CS_trace{b},bin_size,2); CS=CS(:,[1:bin_size:end])./totalSpike{b};
%plot([1:bin_size:size(FineSpikes,2)]/1000,CS,'color',cmap(b-1,:))
%errorbar([1:bin_size:size(FineSpikes,2)]/1000,mean(CS,1),std(CS,0,1)/sqrt(size(CS,1)),'marker','+','color',cmap(b-1,:))
plot([1:bin_size:size(FineSpikes,2)]/1000,mean(CS,1,'omitnan'),'marker','*','color',cmap(b-1,:))
hold all
xlabel('Time (sec)')
ylabel('Fraction of spike')

ax2=subplot(1,4,2);
title('Single spikes')
SS=movsum(SS_trace{b},bin_size,2); SS=SS(:,[1:bin_size:end]);
%plot([1:bin_size:size(FineSpikes,2)]/1000,SS,'color',cmap(b-1,:))
%errorbar([1:bin_size:size(FineSpikes,2)]/1000,mean(SS,1),std(SS,0,1)/sqrt(size(SS,1)),'marker','+','color',cmap(b-1,:))
plot([1:bin_size:size(FineSpikes,2)]/1000,mean(SS,1,'omitnan'),'marker','*','color',cmap(b-1,:))
hold all
xlabel('Time (sec)')
ylabel('Number of spike')

ax4=subplot(1,4,4);
title('Single spikes')
SS=movsum(SS_trace{b},bin_size,2); SS=SS(:,[1:bin_size:end])./totalSpike{b};
%plot([1:bin_size:size(FineSpikes,2)]/1000,SS,'color',cmap(b-1,:))
%errorbar([1:bin_size:size(FineSpikes,2)]/1000,mean(SS,1),std(SS,0,1)/sqrt(size(SS,1)),'marker','+','color',cmap(b-1,:))
plot([1:bin_size:size(FineSpikes,2)]/1000,mean(SS,1,'omitnan'),'marker','*','color',cmap(b-1,:))
hold all
xlabel('Time (sec)')
ylabel('Fraction of spike')
%CS=movsum(CS_trace{b},bin_size,2); CS=CS(:,[1:bin_size:end]);%./sum(FineSpikes,2);
%SS=movsum(SS_trace{b},bin_size,2); SS=SS(:,[1:bin_size:end]);%./sum(FineSpikes,2);
%N=size(FineSpikes,1);
%errorbar([1:bin_size:size(Transient_traces,2)]/1000,mean(CS,1),std(CS,0,1)/sqrt(N),'marker','+','color',cmap(1,:))
end
end
linkaxes([ax1 ax2],'xy')
linkaxes([ax3 ax4],'xy')

%ylabel('Fraction of spike')

%%
AvailResult=1-cellfun(@isempty,Result);
cmap=turbo(size(NiceNeuronList,2));
clear Adap_data
figure;
for i=1:size(NiceNeuronList,2) % Cellc 
    gg=1;
    Av=find(AvailResult(NiceNeuronList(1,i),:));
    for j=[NiceNeuronList(2,i) setdiff(Av,NiceNeuronList(2,i))] % Sessions
        bw=bwlabel(Result{NiceNeuronList(1,i),j}.blue>0);
        if max(bw)>3 
            delay=20; g=1; b=1;
            while g
                blue_seg=find(bw==b);
                if sum(Result{NiceNeuronList(1,i),j}.spike(blue_seg(1):blue_seg(end)+delay))>0
                    g=0;
                    Adap_data{i,1}=Result{NiceNeuronList(1,i),j}.blue(blue_seg(1));
                end
                b=b+1;
            end

            for b=max(bw)-5:max(bw)-1
                blue_seg=find(bw==b);
                Adap_data{i,2}(gg,1)=Result{NiceNeuronList(1,i),j}.blue(blue_seg(1));
                Nss=sum(Result{NiceNeuronList(1,i),j}.spike(blue_seg(1):blue_seg(1)+500));
                Ncs=sum(Result{NiceNeuronList(1,i),j}.Complex(blue_seg(1):blue_seg(1)+500) & Result{NiceNeuronList(1,i),j}.spike(blue_seg(1):blue_seg(1)+500));
                Adap_data{i,2}(gg,2:3)=[Nss Ncs];
                gg=gg+1;
            end

        end
    end % Sessions
Adap_data{i,3}=Adap_data{i,2}; 
denom=   Blue_CAL(find_index_bh(round(Blue_CAL(:,1)*100)/100,round(Adap_data{i,1}(1,1)*100)/100),2);
Norm_int=Blue_CAL(find_index_bh(round(Blue_CAL(:,1)*100)/100,round(Adap_data{i,2}(:,1)*100)/100),2)/denom;
Adap_data{i,3}(:,1)=Norm_int;
plot(Adap_data{i,3}(:,1),Adap_data{i,3}(:,3)./Adap_data{i,3}(:,2),'.','color',cmap(2,:),'markersize',18)
hold all
end % Cell

dat=cell2mat(Adap_data(:,3));




%%

for i=1%:length(fpath)
    load([fpath{i} '/output_data.mat'])
    sz=Device_Data{1, 4}.ROI([2 4]);

    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
    mov_res= mov_mc-mean(mov_mc,3);
    im_G=imgaussfilt(mean(mov_mc,3),3);
    [centers radii]=Cell_segment_circle_10x(im_G,0.7);
    Result{i}.centers=cell_detection_manual(mean(mov_mc,3),centers,[]);
    Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],20);
    Result{i}.traces=-(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt))';

end



%%
for i=1:length(fpath)
    show_traces_spikes(Result{i}.traces,Result{i}.spike,Result{i}.spike)
end
%%
mov_mc=double(mov_mc);
[roi sig]=clicky(mov_mc);
mov_res=mov_mc-movmean(mov_mc,300,3);
%%
%sig=-sig';
tmp=(int) - movmedian(int,70,2); tmp=tmp./get_threshold(tmp,1);
spike=find_spike_bh(tmp,3.5);
plot(int)
hold all
S=int; S(~spike)=NaN;
plot(S,'r.')

%%
ss=find(spike);
for s=1:length(ss)
    try
        sta_mov(:,:,:,s)=mov_res(:,:,ss(s)-50:ss(s)+100);
    end
end
STA=mean(sta_mov,4);
%%
for i=1:length(fpath)
    load([fpath{i} '/output_data.mat'])
    sz=Device_Data{1, 4}.ROI([2 4]);
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz1,sz2));
end
%% High-pass filter
mov_mc=double(mov);
im_G=imgaussfilt(mean(mov_mc,3),2);
mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[8 8]),2);
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[16 16]),2);
[centers radii]=Cell_segment_circle_10x(im_G);
centers=cell_detection_manual(mean(mov_mc,3),centers);
%%
mcTrace = squeeze(mean(xyField,[1 2]));
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));
c_ftprnt=mask_footprint(centers,mov_res,[],10);
show_footprnt(c_ftprnt,mov_mc)
coord=get_coord(c_ftprnt);

traces=-(tovec(mov_mc)'*tovec(c_ftprnt>0))';
traces_hi=squeeze(traces) - movmean(squeeze(traces),400,2);
%traces_hi=squeeze(traces) - mean(squeeze(traces),2);
traces_hi=traces_hi./range(traces_hi,2);
%traces_hi=-traces_hi./mean(traces_hi,2);
%%
figure; scale=0.6;
tiledlayout(10,4)
ax1 = nexttile([2 2]);
colr = max(colormap(jet(size(c_ftprnt,3))),0);
imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

ax4 = nexttile([2 2]);
imshow2(mean(mov_mc,3),[]); colormap('gray')

linkaxes([ax1 ax4],'xy')

ax2 = nexttile([7 4]);
lines=plot(traces_hi'+[1:size(traces_hi,1)]*scale);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(size(c_ftprnt,3)),2))
axis tight
set(gca,'ytick',[1:size(traces_hi,1)]*scale,'yticklabel',[1:size(traces_hi,1)])
ax3 = nexttile([1 4]);
plot(Blue)

linkaxes([ax2 ax3],'x')
saveas(gcf,[fpath 'voltage_trace_plot'])
save([fpath 'result.mat'],'c_ftprnt','traces')

%% Complex spike -> single spike adaptation figure generation

% load data
clear
cd('/Users/bhlee1117/Dropbox/20230718')
save_to='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Manuscript/2023_NaV_InactivationFigure';
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

V2Result(8,:)=V2Result(8,[1 3 2]);
V2Result(2,:)=V2Result(2,[1 3 2]);

V2Result([2 3 15],:)=[];
%% find complex spike and transients
clear Transients 

for i=1:size(V2Result,1) %neurons
    for Sessions=1:size(V2Result,2)
        if ~isempty(V2Result{i,Sessions})
            tr=V2Result{i,Sessions}.trace;
            blue=bwlabel(V2Result{i,Sessions}.Blue);
          
            tr=tr-median(tr); tr=tr./get_threshold(tr,1); %normalize to noise level
            tr_CS=tr-movmedian(tr,300,2); % high pass
            tr_CS=tr_CS/get_threshold(tr_CS,1); % normalize to noise level
            tr_CS2=tr-median(tr,2); % high pass
            for b=1:max(blue)
                if sum(blue==b)>400 & sum(blue==b)<600
                tr_CS2(find(blue==b))=tr_CS2(find(blue==b))-median(tr_CS2(find(blue==b)));
                else
                tr_CS2(find(blue==b))=tr_CS(find(blue==b));    
                end
            end

            tr_SS=tr-movmedian(tr,30,2); % high pass
            tr_SS=tr_SS/get_threshold(tr_SS,1); % normalize to noise level

            spike_SS=find_spike_bh(tr_SS,4,3);
            spike_complex=find_spike_bh(tr_CS,4,2);

            [trans tr_trace]=detect_transient(movmean(tr_CS2,5,2),[3 1]); %transient detection
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
            spike_total=zeros(1,size(tr,2));
            spike_total(find(spike_SS))=1;
            spike_total(find(spike_complex))=2;
            V2Result{i,Sessions}.spike=spike_total;
            V2Result{i,Sessions}.CSTransient=ComplexTransient>0;

        end
    end
end


%% Show representative Trace of pulse stimulation
noi=3; DAQdT=1e-5;

load([V2Result{noi,1}.fpath '/output_data.mat'])
camTrig=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
t=camTrig*DAQdT;

f1=figure;
tiledlayout(3,1)
ax1=nexttile([2 1]);
CSTr=V2Result{noi,1}.trace; CSTr(V2Result{noi,1}.CSTransient==0)=NaN;
plot(t,V2Result{noi,1}.trace,'linewidth',1,'color',[0 0.1 0.3])
hold all
plot(t,CSTr,'linewidth',1,'color',[0.7 0 0.3])
sp=find(V2Result{noi,1}.spike);
plot(t(sp),V2Result{noi,1}.trace(sp),'r.','markersize',8)
axis off

ax2=nexttile([1 1]);
plot(t,V2Result{noi,1}.Blue,'linewidth',3)
axis off

linkaxes([ax1 ax2],'x')
xlim([3.5 9])

set(f1, 'Position', [100, 100, 1000, 300]);

%% Show representative Trace of Ramp stimulation
noi=3; DAQdT=1e-5;

load([V2Result{noi,1}.fpath '/output_data.mat'])
camTrig=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
t=camTrig*DAQdT;

f1=figure;
tiledlayout(3,1)
ax1=nexttile([2 1]);
CSTr=V2Result{noi,2}.trace; CSTr(V2Result{noi,2}.CSTransient==0)=NaN;
plot(t,V2Result{noi,2}.trace,'linewidth',1,'color',[0 0.1 0.3])
hold all
plot(t,CSTr,'linewidth',1,'color',[0.7 0 0.3])
sp=find(V2Result{noi,2}.spike);
plot(t(sp),V2Result{noi,2}.trace(sp),'r.','markersize',8)
axis off

ax2=nexttile([1 1]);
plot(t,V2Result{noi,2}.Blue,'linewidth',3)
axis off

linkaxes([ax1 ax2],'x')
xlim([0.5 6.5])

set(f1, 'Position', [100, 100, 1000, 300]);

%% Show Statistics (pulse)
session=1; DAQdT=1e-5;
FracCS=zeros(size(V2Result,1),5);
FracCS_bin=zeros(size(V2Result,1),5,4);
for i=1:size(V2Result,1)
    t_segment=[3.5 9];
    load([V2Result{i,session}.fpath '/output_data.mat'])
    camTrig=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    t=camTrig*DAQdT;
    toi=find(t>t_segment(1) & t<t_segment(2));
            
    tr=V2Result{i,session}.trace(toi);
    sp=V2Result{i,session}.spike(toi);
    bl=V2Result{i,session}.Blue(toi);               

    LabelBlue=bwlabel(bl);
    
    for l=1:5
        Blue_index=find(LabelBlue==l);
        NumberCS=sum(sp(Blue_index)==2);
        NumberSpike=sum(sp(Blue_index)>0);
        FracCS(i,l)=NumberCS/NumberSpike;
        
        pulse_bin=[1:floor(length(Blue_index)/4):length(Blue_index)];
        for p=1:length(pulse_bin)-1
        FracCS_bin(i,l,p)=sum(sp(Blue_index(pulse_bin(p):pulse_bin(p+1)-1))==2)/sum(sp(Blue_index(pulse_bin(p):pulse_bin(p+1)-1))>0);
        end
    end
end
    f1=figure;      
    plot([1:5],FracCS','color',[0.6 0.6 0.6],'marker','.','linewidth',1)
    hold all
    errorbar([1:5],mean(FracCS,1,'omitnan'),std(FracCS,0,1,'omitnan')/sqrt(size(V2Result,1)),'color',[0.7 0 0.3],'marker','+','linewidth',3)
    xlim([0.5 5.5])
    xlabel('Blue light strength')
    ylabel('Fraction of Complex Spikes')
    set(gca,'Fontsize',15)
    set(f1, 'Position', [100, 100, 300, 300]);
    saveas(f1,[save_to '/Pulse_CSratio_Blue.fig'])
    print(f1, [save_to '/Pulse_CSratio_Blue.jpg'],'-djpeg', ['-r', num2str(600)]);
    
    f1=figure;      
    cmap=jet(5);
    for p=1:5
    errorbar([1:4],squeeze(mean(FracCS_bin(:,p,:),1,'omitnan')),squeeze(std(FracCS_bin(:,p,:),0,1,'omitnan'))/sqrt(size(V2Result,1)), ...
        'color',cmap(p,:),'marker','+','linewidth',3)
    hold all
    end
    xlim([0.5 4.5])
    xlabel('Time (ms)')
    ylabel('Fraction of Complex Spikes')
    set(gca,'xtick',[1:4],'xticklabel',[1:4]*125,'Fontsize',15)
    set(f1, 'Position', [100, 100, 300, 300]);
    saveas(f1,[save_to '/Pulse_CSratio_Time.fig'])
    print(f1, [save_to '/Pulse_CSratio_Time.jpg'],'-djpeg', ['-r', num2str(600)]);

%% Show Statistics (Ramp)
session=2; DAQdT=1e-5;
Bin_number=30;    
CS_bin=zeros(size(V2Result,1),Bin_number); Spike_bin=zeros(size(V2Result,1),Bin_number);
load([V2Result{5,1}.fpath '/output_data.mat'])
camTrig=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
t=camTrig*DAQdT;
for i=1:size(V2Result,1)
    
    t_segment=t(bwlabel(V2Result{i,session}.Blue)==1);
    
    toi=find(t>t_segment(1) & t<t_segment(end));
            
    tr=V2Result{i,session}.trace(toi);
    sp=V2Result{i,session}.spike(toi);
    bl=V2Result{i,session}.Blue(toi);               

    LabelBlue=bwlabel(bl);
    
        Blue_index=find(LabelBlue==1);
        NumberCS=sum(sp(Blue_index)==2);
        NumberSpike=sum(sp(Blue_index)>0);
        
        pulse_bin=[1:floor(length(Blue_index)/Bin_number):length(Blue_index)];
        for p=1:length(pulse_bin)-1
        CS_bin(i,p)=sum(sp(Blue_index(pulse_bin(p):pulse_bin(p+1)-1))==2);
        Spike_bin(i,p)=sum(sp(Blue_index(pulse_bin(p):pulse_bin(p+1)-1))>0);
        end
end
    
        f1=figure;      
        tiledlayout(2,1)
        ax1=nexttile([1 1]);
        bar([1:Bin_number],sum(Spike_bin,1),'FaceColor',[0.4 0.4 0.4]);
        hold all
        bar([1:Bin_number],sum(CS_bin,1),'FaceColor',[0.7 0 0.3]);
        ylabel('Number of Spikes')
        set(gca,'xtick',[0:5:Bin_number],'xticklabel',[0:5:Bin_number]*6000/Bin_number,'Fontsize',15)
    
        ax2=nexttile([1 1]);
        errorbar([1:Bin_number],mean(CS_bin./Spike_bin,1,'omitnan'),std(CS_bin./Spike_bin,0,1,'omitnan')/sqrt(size(V2Result,1)), ...
            'color',[0.7 0 0.3],'marker','+','linewidth',3);
        hold all
        errorbar([1:Bin_number],mean((Spike_bin-CS_bin)./Spike_bin,1,'omitnan'),std((Spike_bin-CS_bin)./Spike_bin,0,1,'omitnan')/sqrt(size(V2Result,1)), ...
            'color',[0 0.1 0.3],'marker','+','linewidth',3);
        linkaxes([ax1 ax2],'x')
        xlim([0.5 31])
        
        xlabel('Time after stimulation onset (ms)')
        ylabel('Fraction of Complex Spikes')
        set(gca,'xtick',[0:5:Bin_number],'xticklabel',[0:5:Bin_number]*6000/Bin_number,'Fontsize',15)
        set(f1, 'Position', [100, 100, 600, 300]);
        saveas(f1,[save_to '/Ramp_CSratio_Time.fig'])
        print(f1, [save_to '/Ramp_CSratio_Time.jpg'],'-djpeg', ['-r', num2str(600)]);

%% sort complex spikes


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
% Wavelet transform
Fs = 1000; % Sampling frequency (change as per your signal)
t = 0:1/Fs:1; % Time vector

waveletName = 'cmor'; % Wavelet family (you can choose other wavelets as well)
scales = 1:128; % Range of scales (adjust as needed)

[coefficients{i}, frq{i}] = cwt(tr, 1000);
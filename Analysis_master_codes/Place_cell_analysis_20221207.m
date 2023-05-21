
%addpath(genpath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code'));
clear
[fpath] = uigetfile_n_dir();

%%
for i=4%3:length(fpath)

    Sz = importdata([fpath{i} '/experimental_parameters.txt']);
    load(fullfile(fpath{i},'settings.mat'));
    sz1=Sz.data(1); sz2=Sz.data(2);
    [a,b]=system(sprintf('GetFileInfo "%s"',fullfile(fpath{i},'settings.mat'))); s=strfind(b,'modified:')+10; crdat=b(s:s+18);
    image_start=datestr(datenum(crdat));
    dt=datetime(image_start(end-8:end),'InputFormat','HH:mm:ss'); dt.Format = 'HH:mm:ss';
    fnm=dir(fullfile(fpath{i}, '*.csv'));

    a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);
    frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];
    [Arduino_data reward_pos lap_dist]=match_treadmill_DAQ(fullfile(fpath{i},fnm(2).name),1e-5,DAQ_data,(1/frm_rate+0.00002),2); % time, treadmill, Reward, Run
    Arduino_data(:,1)=Arduino_data(:,1)*(1/frm_rate*1000)/(1/frm_rate*1000+0.02);
    aa=find(bwlabel(Arduino_data(:,4))==1);

    frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];

    mov_mc=double(readBinMov_times([fpath{i} '/mc01.bin'],sz2,sz1,[aa(1):aa(1)+500]));
    im_G=imgaussfilt(mean(mov_mc,3),2);
    [centers radii]=Cell_segment_circle_2x(im_G,0.9);
     %Result{i}.centers=cell_detection_manual(mean(mov_mc,3),centers,[0 9000]);
    Result{i}.centers=cell_detection_manual(mean(mov_mc,3),Result{3}.centers,[0 9000]);

end

%%
for i=4%length(fpath)
    load(fullfile(fpath{i},'settings.mat'));
    load(fullfile(fpath{i},'mcTrace01.mat'));
    Sz = importdata([fpath{i} '/experimental_parameters.txt']);
    sz1=Sz.data(1); sz2=Sz.data(2);

    [a,b]=system(sprintf('GetFileInfo "%s"',fullfile(fpath{i},'settings.mat'))); s=strfind(b,'modified:')+10; crdat=b(s:s+18);
    image_start=datestr(datenum(crdat));
    dt=datetime(image_start(end-8:end),'InputFormat','HH:mm:ss'); dt.Format = 'HH:mm:ss';
    fnm=dir(fullfile(fpath{i}, '*.csv'));

    a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);
    frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];
    [Arduino_data reward_pos lap_dist]=match_treadmill_DAQ(fullfile(fpath{i},fnm(2).name),1e-5,DAQ_data,(1/frm_rate+0.00002),2); % time, treadmill, Reward, Run
    Arduino_data(:,1)=Arduino_data(:,1)*(1/frm_rate*1000)/(1/frm_rate*1000+0.02);
    Result{i}.Arduino=Arduino_data;

    mov_mc=double(readBinMov([fpath{i} '/mc' num2str(1,'%02d') '.bin'],sz2,sz1));
    mov_res= mov_mc(:,:,1:10000)-mean(mov_mc,3);

    mov_res = SeeResiduals(mov_res,mcTrace(1:10000,:));
    mov_res = SeeResiduals(mov_res,mcTrace(1:10000,:).^2);
    mov_res = SeeResiduals(mov_res,mcTrace(1:10000,1).*mcTrace(1:10000,2));

    Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6);
    Result{i}.ref_im=mean(mov_mc,3);
    Result{i}.traces=[]; Result{i}.traces_hi=[]; Result{i}.traces_bin=[]; Result{i}.traces_bin_hi=[];
    Result{i}.traces_res=[]; Result{i}.traces_res_hi=[]; Result{i}.mcTrace=[];
    Result{i}.im_corr=[];

    Result{i}.frm_rate=frm_rate;
    Result{i}.DMDPatt=dmd_mask_sequence_rois;
    Result{i}.coord=get_coord(Result{i}.c_ftprnt);
    Result{i}.traces=[Result{i}.traces -(tovec(mov_mc(:,:,1:10000))'*tovec(Result{i}.c_ftprnt))'];
    Result{i}.traces_bin=[Result{i}.traces_bin -(tovec(mov_mc(:,:,1:10000))'*tovec(Result{i}.c_ftprnt>0))'];
    mov_mc_vec=tovec(mov_mc(:,:,1:10000)); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ref_im_vec=tovec(Result{i}.ref_im); ref_im_vec=(ref_im_vec-mean(ref_im_vec))./std(ref_im_vec);
    Result{i}.im_corr=sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1);
    %Result{i}.traces_res=[Result{i}.traces_res -(tovec(mov_res(:,:,1:10000))'*tovec(Result{i}.c_ftprnt))'];

    for j=2:length(f_seg)-1
        mov_mc=double(readBinMov([fpath{i} '/mc' num2str(j,'%02d') '.bin'],sz2,sz1));
        try mov_mc=mov_mc(:,:,1:10000); catch mov_mc=mov_mc; end

        load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);
        %mov_res= mov_mc-mean(mov_mc,3);
        %mov_res = SeeResiduals(mov_res,mcTrace(1:size(mov_mc,3),:));
        %mov_res = SeeResiduals(mov_res,mcTrace(1:size(mov_mc,3),:).^2);
        %mov_res = SeeResiduals(mov_res,mcTrace(1:size(mov_mc,3),1).*mcTrace(1:size(mov_mc,3),2));
        mov_mc_vec=tovec(mov_mc); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
        Result{i}.im_corr=[Result{i}.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

        Result{i}.traces=[Result{i}.traces -(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt))'];
        %Result{i}.traces_res=[Result{i}.traces_res -(tovec(mov_res)'*tovec(Result{i}.c_ftprnt))'];
        Result{i}.traces_bin=[Result{i}.traces_bin -(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt>0))'];
        j
    end
    for j=1:length(f_seg)-1
        load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);
        try mcTrace=mcTrace(1:10000,:); catch mcTrace=mcTrace; end
        Result{i}.mcTrace=[Result{i}.mcTrace; mcTrace];
    end

    [V, D, W] = eig(corrcoef(zscore(Result{i}.traces,0,2)'));
    D = diag(D); D = D(end:-1:1); V = V(:,end:-1:1);
    PCA_trace=zscore(V(:,1)'*Result{i}.traces,0,2);
    ns=size(Result{i}.traces,1);

    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces,ns,1,[]),Result{i}.mcTrace'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),Result{i}.mcTrace.^2'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),(Result{i}.mcTrace(:,1).*Result{i}.mcTrace(:,2))'));
    %Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),PCA_trace'));
    Result{i}.traces_res=moving_residual(Result{i}.traces_res,Result{i}.mcTrace,PCA_trace,300);
    Result{i}.traces_res_hi=Result{i}.traces_res - movmedian(Result{i}.traces_res,150,2);

    Result{i}.traces_hi=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),150,2);
    Result{i}.traces_bin_hi=squeeze(Result{i}.traces_bin) - movmedian(squeeze(Result{i}.traces_bin),200,2);

    Result{i}.noise=get_threshold(Result{i}.traces,1); Result{i}.noise_res=get_threshold(Result{i}.traces_res,1);
    Result{i}.traces_hi=Result{i}.traces_hi./Result{i}.noise;
    Result{i}.traces_bin_hi=Result{i}.traces_bin_hi./get_threshold(Result{i}.traces_bin_hi,1);
    Result{i}.traces_res_hi=Result{i}.traces_res_hi./Result{i}.noise_res;

    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),20,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=find_spike_bh(tmp,3,1.5);
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),100,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3,1.5))>0;
end
%save(['Treadmill7-9_20230209_result.mat'],'Result','fpath','-v7.3')
save(['Treadmill_BHLm008_20230209_result.mat'],'Result','fpath','-v7.3')
%%
% n=81;
%
% S_trace=zscore(Result{i}.traces,0,2);
% covMat=corrcoef(S_trace');
% [V, D, W] = eig(covMat);
%
%     D = diag(D);
%     D = D(end:-1:1);
%     V = V(:,end:-1:1);
%     vSign = sign(max(V) - max(-V));  % make the largest value always positive
%     V = V.*vSign;
% pca_trace=V(:,1)'*S_trace;
% res_trace=squeeze(SeeResiduals(reshape(Result{i}.traces,135,1,[]),Result{i}.mcTrace'));
% res_trace=squeeze(SeeResiduals(reshape(res_trace,135,1,[]),Result{i}.mcTrace.^2'));
% res_trace=squeeze(SeeResiduals(reshape(res_trace,135,1,[]),(Result{i}.mcTrace(:,1).*Result{i}.mcTrace(:,2))'));
% res_trace2=squeeze(SeeResiduals(reshape(res_trace,135,1,[]),pca_trace));
% %res_trace3=squeeze(SeeResiduals(reshape(res_trace,135,1,[]),V(:,1:2)'*S_trace));
% res_trace3=moving_residual(Result{i}.traces_res,Result{i}.mcTrace,pca_trace);
% tiledlayout(9,1)
% ax1=nexttile([7 1]);
% plot(Result{i}.traces(n,:)-movmedian(Result{i}.traces(n,:),200,2))
% hold all
% plot(res_trace(n,:)-movmedian(res_trace(n,:),200,2)+13000)
% plot(res_trace2(n,:)-movmedian(res_trace2(n,:),200,2)+26000)
% plot(res_trace3(n,:)-movmedian(res_trace3(n,:),200,2)+39000)
% axis tight
% %plot(Result{i}.traces_res(n,:)-movmedian(Result{i}.traces_res(n,:),200,2)+26000)
% ax2=nexttile([2 1]);
% plot(Result{i}.mcTrace)
% yyaxis right
% plot(zscore(V(:,1)'*S_trace,0,2))
% axis tight
% linkaxes([ax1 ax2],'x')
%%


%%
figure; i=4;
tiledlayout(11,1)
ax1=nexttile([7 1]);
noi=[1:size(Result{i}.traces_hi,1)];
%noi=[1:50];
t=[1:size(Result{i}.traces_hi,2)]/800; scale=10;
%tr=Result{i}.traces_res-median(Result{i}.traces_res,2); fprnt=Result{i}.c_ftprnt;
tr=Result{i}.traces_res-movmedian(Result{i}.traces_res,200,2);
tr=tr./get_threshold(tr,1);

lines=plot(t,tr(noi,:)'+[1:size(noi,2)]*scale); line_color=jet(length(noi));
%arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
hold all
S=tr; S(~Result{i}.spike)=NaN;
plot(t,S(noi,:)'+[1:size(noi,2)]*scale,'r.')
set(gca,'ytick',[1:size(noi,2)]*scale,'yticklabel',noi)

ax2=nexttile([2 1]);
plot(Result{i}.Arduino(:,1),Result{i}.Arduino(:,2))
hold all
d=Result{i}.Arduino(:,2); d(~Result{i}.Arduino(:,4))=NaN;
plot(Result{i}.Arduino(:,1),d)
plot(Result{i}.Arduino(:,1),Result{i}.Arduino(:,3)*5000,'r')


ax3=nexttile([2 1]);
plot(Result{i}.Arduino(:,1),Result{i}.mcTrace','k')
hold all
d=Result{i}.mcTrace; d(~Result{i}.Arduino(:,4),:)=NaN;
plot(Result{i}.Arduino(:,1),d','r')
axis tight

linkaxes([ax1 ax2 ax3],'x')
%%



%% place voltage movie
noi=[4 6 12 24 26 30 37 51 55];
place_cell_movie(['test_pc_video'],Arduino_data,volt,spike,noi,15000,mov_test,ref_mask)

%% Place triggered average movie
place_bin=30; %number of bin, 45 -> ~20 mm per bin
bin_dist=ceil(Arduino_data(:,2)/(lap_dist/place_bin)); bin_dist(~Arduino_data(:,4))=NaN;
for p=1:place_bin; PTA{p}=[]; end

for i=1:length(f_seg)-1
    i
    % load movie during walking
    mov_seg=double(readBinMov_times([fpath '/motion_corrected/mc' num2str(i,'%02d') '.bin']...
        ,size(mov_test,1),size(mov_test,2),find(Arduino_data(f_seg(i):f_seg(i+1)-1,4))));
    tmp=bin_dist(f_seg(i):f_seg(i+1)-1); tmp=tmp(find(Arduino_data(f_seg(i):f_seg(i+1)-1,4)));
    p_l=(unique(bin_dist(f_seg(i):f_seg(i+1)-1)))';
    for p=p_l(find(~isnan(p_l)))
        PTA{p}(:,:,[size(PTA{p},3)+1:size(PTA{p},3)+sum(tmp==p)])=mean(mov_seg,3)-(mov_seg(:,:,find(tmp==p)));
    end
end
for p=1:place_bin
    PTA_movie(:,:,p)=mean(PTA{p},3,'omitnan');
end

%% Place triggered average trace; Voltage trace, Spike
f=1;
place_bin=30; %number of bin, 45 -> ~20 mm per bin
bin_dist=ceil(Result{f}.Arduino(:,2)/(mean(lap_dist)/place_bin));
bin_dist(~Result{f}.Arduino(:,4))=NaN;
clear place_volt place_volt_lp place_spike
for i=1:size(Result{f}.traces_hi,1)
    for p=1:place_bin
        tmp=find(bin_dist==p);
        place_volt(i,p)=mean(Result{f}.traces(i,tmp));
        place_volt_lp(i,p)=mean(Result{f}.traces_hi(i,tmp));
        place_spike(i,p)=mean(Result{f}.spike(i,tmp));
    end
    %place_volt=movmean(place_volt,3,2);
    place_volt=(place_volt-min(place_volt,[],2))./(max(place_volt,[],2)-min(place_volt,[],2));
    %place_volt_lp=movmean(place_volt_lp,3,2);
    place_volt_lp=(place_volt_lp-min(place_volt_lp,[],2))./(max(place_volt_lp,[],2)-min(place_volt_lp,[],2));
    place_spike=movmean(place_spike,3,2);%-min(place_volt,2);
    place_spike=(place_spike-min(place_spike,[],2))./(max(place_spike,[],2)-min(place_spike,[],2));
end

% sort by spike #
[a b]=sort(place_spike,2,'descend');
put_pc=find(max(place_spike,[],2)>std(a(:,1:round(size(a,2)*0.3)),0,2)*2);
%calculate centroid of place field
%[a c]=sort(place_spike(put_pc,:)*[1:place_bin]'./sum(place_spike(put_pc,:),2),'ascend'); %calculate centroid
%calculate maximum of place field
[a arg]=max(place_spike(put_pc,:),[],2);
[a c]=sort(arg);
order=[c ;setdiff([1:size(Result{f}.traces_hi,1)],put_pc)'];
figure;
subplot(1,3,1)
imagesc(place_volt(order,:)); axis tight; hold on;
line(repmat(ceil(reward_pos/(mean(lap_dist)/place_bin)),1,2),[0 length(order)+0.5],'color','r');
subplot(1,3,2)
imagesc(place_volt_lp(order,:)); axis tight; hold on;
line(repmat(ceil(reward_pos/(mean(lap_dist)/place_bin)),1,2),[0 length(order)+0.5],'color','r');
subplot(1,3,3)
imagesc(place_spike(order,:)); axis tight; hold on;
line(repmat(ceil(reward_pos/(mean(lap_dist)/place_bin)),1,2),[0 length(order)+0.5],'color','r');

%% now plot by every laps
clear lap_place_spike lap_place_volt lap_place_volt_lp
figure; noi=[1:10]; l=bwlabel(ceil(Arduino_data(:,2)/(mean(lap_dist)/place_bin))==place_bin)+1;
g=1; R=zeros(size(l,1),1);
show_row=2;
for i=1:max(l); [a b]=max(find(l==i)); R(a)=1; end; R=cumsum(R);
for n=noi
    for r=1:max(R)-1 %laps
        for p=1:place_bin %place bin
            tmp=(bin_dist==p & R==r);
            if length(find(tmp))==1; tmp=[]; end
            lap_place_volt{n}(r,p)=mean(Result{f}.traces_hi(n,find(tmp)),'omitnan');
            lap_place_spike{n}(r,p)=mean(Result{f}.spike(n,find(tmp)),'omitnan');
        end
    end
    % voltage trace
    %subplot(show_row,16,(ceil(g/8)-1)*16+mod(g,8)*2-1+16*(mod(g,8)==0))
    nexttile;
    %figure;
    imagesc([lap_place_volt{n}-min(lap_place_volt{n}(:)); zeros(3,place_bin);...
        10*place_volt(n,:)]);

    hold on;
    line(repmat(ceil(reward_pos/(mean(lap_dist)/place_bin)),1,2),[0 max(R)+0.5],'color','r','linestyle',':');
    text(1,-5,[num2str(n) 'th neuron'],'Fontsize',6)
    axis tight off
    %figure;
    %plot([1:place_bin],place_volt(n,:))
    % mean voltage trace by place bin
    %figure;
    nexttile;
    imagesc([lap_place_spike{n}; zeros(3,place_bin);...
        0.1*(place_spike(n,:)-min(place_spike(n,:)))/(max(place_spike(n,:)-min(place_spike(n,:))))]);
    hold on;
    line(repmat(ceil(reward_pos/(mean(lap_dist)/place_bin)),1,2),[0 max(R)+0.5],'color','r','linestyle',':');
    axis tight off
    %figure;
    %plot([1:place_bin],place_spike(n,:))
    g=g+1;
end


%% show raw voltage trace at place bin
cmap=winter(20); place_bin=30;
bin_dist=ceil(Result{f}.Arduino(:,2)/(mean(lap_dist)/place_bin));
l=bwlabel(ceil(Result{f}.Arduino(:,2)/(mean(lap_dist)/place_bin))==place_bin)+1; R=zeros(size(l,1),1);
for i=1:max(l); [a b]=max(find(l==i)); R(a)=1; end; R=cumsum(R);
figure; noi=75; rng=[5000 3000]; scale=10;
for i=0:max(R)
    try
        tmp=find((bin_dist==29 & R==i));
        ind(i+1)=tmp(1);
        plot([1:sum(rng)+1]*1.25,Result{f}.traces_hi(noi,tmp(1)-rng(1):tmp(1)+rng(2))+i*scale,'color',cmap(i+1,:))
        hold all
        line([1 sum(rng)+1]*1.25,[i*scale i*scale],'color',[0.6 0.6 0.6]);
        text([-400]*1.25,i*scale,[num2str((tmp(1)-rng(1))*1.25*1e-3,'%.1f') 's'])
        text([sum(rng)+50]*1.25,i*scale,[num2str((tmp(1)+rng(2))*1.25*1e-3,'%.1f') 's'])
        s_tmp=find(Result{f}.spike(noi,tmp(1)-rng(1):tmp(1)+rng(2)));
        plot([s_tmp]*1.25,Result{f}.traces_hi(noi,tmp(1)-rng(1)-1+s_tmp)+i*scale,'r.')
    end
end
yyaxis right
line([rng(1) rng(1)]*1.25,[0 1],'color','c')
axis off
% gg=1;
% for p=[25:30]
%     t=nan(3001,1); tmp2=bin_dist==p;
%     plot([1:1:4001]*1.25,tmp2(tmp(1)-3000:tmp(1)+1000),'color',cmap(gg,:),'marker','.')
%     gg=gg+1;
% end
% axis off



%% Convert to firing rate and plot Firing rate vs place

bin_time=500; %ms
tmp=movsum(spike,bin_time/frm_rate,2);
frame_ind=[bin_time/frm_rate-1:bin_time/frm_rate:size(tmp,2)];
%FR=tmp(:,frame_ind)*1000/bin_time; %Firing rate (Hz)
FR=tmp*1000/bin_time;
bin_dist=ceil(Arduino_data(:,2)/(lap_dist/45));

bin_time=100; %different bin time with firing rate
frame_ind=[bin_time/frm_rate-1:bin_time/frm_rate:size(tmp,2)];
FR=FR(:,frame_ind); FR(:,frame_ind<16000)=0;
neuron_id=[38]; g=1;
for i=neuron_id
    cline([bin_time:bin_time:size(FR,2)*bin_time]/1000,g*5000+bin_dist(frame_ind)',...
        zeros(size(FR,2),1)',round(FR(i,:))+1,jet(max(round(FR(i,:))+1)),3);
    %cline([1:1:size(FR,2)],g*5000+Arduino_data(:,2)',...
    %zeros(size(FR,2),1)',round(FR(i,:))+1,jet(max(round(FR(i,:))+1)));
    hold all
    g=g+1;
end
yyaxis right
plot(Arduino_data(:,1),Arduino_data(:,3),'--')
axis tight off
%% Reward/Stim triggered movie
clear rng RTM
range=[100 100];
%R=find(Arduino_data(:,3)==1); g=1;
R=find([zeros(1,20000) spike(37,20001:end)]==1)'; g=1;
for i=R'
    rng(g,:)=[ceil((i-range(1))/(f_seg(2)-1)) mod(i-range(1),f_seg(2)-1) ...
        ceil((i+range(2))/(f_seg(2)-1)) mod(i+range(2),f_seg(2)-1)];
    if rng(g,1)==rng(g,3)
        tmp=double(readBinMov_times([fpath '/motion_corrected/mc' num2str(rng(g,1),'%02d') '.bin']...
            ,size(mov_test,1),size(mov_test,2),[rng(g,2):rng(g,4)]));
    else
        tmp=double(readBinMov_times([fpath '/motion_corrected/mc' num2str(rng(g,1),'%02d') '.bin']...
            ,size(mov_test,1),size(mov_test,2),[rng(g,2):f_seg(2)-1]));
        tmp2=double(readBinMov_times([fpath '/motion_corrected/mc' num2str(rng(g,3),'%02d') '.bin']...
            ,size(mov_test,1),size(mov_test,2),[1:rng(g,4)]));
        tmp=cat(3,tmp,tmp2);
    end
    RTM(:,:,:,g)=tmp;
    g=g+1;
end
RTM_movie=-mean(RTM,4)+mean(mean(RTM,4),3);

%% STA
clear ST
noi=37; range=[100 100];
%s=find([zeros(1,20000) spike(noi,20001:end)]==1);
s=find([spike(noi,:)]==1);
g=1;
for i=s
    try
        ST(:,:,g)=volt(:,i-range(1):i+range(2));
        g=g+1;
    end
end
STA=mean(ST,3);
for i=37%1:size(volt,1)
    plot([-range(1):1:range(1)]*1.25,STA(i,:))
    hold all
    %text([-100],STA(i,1),num2str(i))
end
%% synchronous spike
clear strace synchro_c_list

spike_threshold=5;
sp=find(sum(spike,1)==spike_threshold);
synchro_c_list=NaN(max(sum(spike,1)),length(sp));
g=1;
for s=1:length(sp)
    synchro_c_list(1:sum(spike(:,sp(s))),s)=find(spike(:,sp(s)));
    try
        strace(:,:,g)=volt(:,sp(s)-200:sp(s)+200)-median(volt(:,sp(s)-200:sp(s)+200),2);
        g=g+1;
    end
end

[synchro_cell]=unique(synchro_c_list);
for i=1:size(synchro_cell,1)
    synchro_cell(i,2)=sum((synchro_c_list(:)==synchro_cell(i)));
end
synchro_cell(isnan(synchro_cell(:,1)),:)=[];

sta=mean(strace,3);
figure;
g=0; scale=500;
for noi=synchro_cell(find(synchro_cell(:,2)<=10),1)'
    plot([-200:200]*1.25,sta(noi,:)+g*scale);
    hold all
    g=g+1;
end



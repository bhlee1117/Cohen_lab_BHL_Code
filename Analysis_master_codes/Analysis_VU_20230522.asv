% Analysis on AAV expression sample and plot, in house YQ201
% 2022/11/09, Byung Hun Lee

clear
[fpath] = uigetfile_n_dir;
%% Motion correction
for i=1:length(fpath)
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 4}.ROI([2 4]));
    mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
    %mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[1:30000]));

    mov_test=mov(:,:,150:250);
    try mov_test = single(mov_test)./single(max(mov_test.data(:)));
    catch disp('change to vm')

        mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));
    [mov_mc,xyField]=optical_flow_motion_correction_LBH(vm(mov),mov_ref,'normcorre');
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
    %mcTrace = squeeze(mean(xyField,[1 2]));
    save([fpath{i} '/mcTrace.mat'],'xyField')
end

%moviefixsc(mov)
%hold all
%plot(dmd_current_roi{1}(:,1),dmd_current_roi{1}(:,2))
%Blue=DAQ_waves.amplitude(4,round([1:size(mov,3)]*1.27*1e-3/1e-5));
%% Signal extraction


C=0; DAQ_rate=0.00001;
for i=78:length(fpath)
    s=split(fpath{i},'/');
    s2=split(s(end),'_'); s2=str2num(s2{2});
    if C==s2
        g=g+1;
    else
        g=1;
    end
    C=s2;

    disp([fpath{i} ' Cell#' num2str(s2) ' ' num2str(g)])
    load([fpath{i} '/output_data.mat'])
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    a=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    frm_rate=(a(2)-a(1))*DAQ_rate;
    sz=Device_Data{1, 4}.ROI([2 4]);
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
    mov_res= mov_mc-mean(mov_mc,3);
    im_G=imgaussfilt(mean(mov_mc,3),3);

    t=[1:size(mov_mc,3)-1]*frm_rate; t(t>length(Blue)*DAQ_rate)=[];
    centers=[]; th=0.1;
    while isempty(centers)
        [centers radii]=Cell_segment_circle_10x(im_G,th);
        th=th+0.01;
    end
    if sum((centers-sz/2).^2)>100
        centers=sz/2;
    end
    Result{s2,g}.fp=fpath{i};
    Result{s2,g}.blue=Blue(round(t/DAQ_rate));
    Result{s2,g}.c_ftprnt=mask_footprint(centers,movmean(mov_res(:,:,1000:end),10,3),[],20);
    Result{s2,g}.traces=-(tovec(mov_mc)'*tovec(Result{s2,g}.c_ftprnt))';
    Result{s2,g}.ref_im=mean(mov_mc,3);
    Result{s2,g}.spike=[];
    tmp=squeeze(Result{s2,g}.traces) - movmedian(squeeze(Result{s2,g}.traces),4,2); tmp=tmp./get_threshold(tmp,1);
    Result{s2,g}.spike=find_spike_bh(tmp,3,3);
    tmp=squeeze(Result{s2,g}.traces) - movmedian(squeeze(Result{s2,g}.traces),100,2); tmp=tmp./get_threshold(tmp,1);
    Result{s2,g}.spike=(Result{s2,g}.spike+find_spike_bh(tmp,3,1.5))>0;
end


%%
clear
load('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230408_BU_BHLm026/20230521_Adaptation_result_BHLm2526.mat')
L=sum(1-cellfun(@isempty,Result),2);
list_pyr=find(L>1)';
%y=[8 16 19 23 27 29 33 34 36 37 41 42 45 46 48 60 61 65 68 69 70 73 83 84 88 90 97 102 104 107 110 114 117 118];

%% Show cell images
figure;
for i=y(1,:)
    nexttile([1 1])
    imshow2(Result{i,1}.ref_im,[])
    nexttile([1 1])
    imshow2(Result{i,1}.c_ftprnt(:,:,1),[])
end

%%
AvailResult=1-cellfun(@isempty,Result);
cmap=winter(size(y,2));
g=1;
figure;
tick=[];

for i=30%1:size(y,2)
    ax1=[];
    %ax1=[ax1 nexttile([1 1])];
    Av=find(AvailResult(y(1,i),:));

    for j=[y(2,i) setdiff(Av,y(2,i))]
        t=[1:length(Result{y(1,i),j}.blue)];
        tr=Result{y(1,i),j}.traces(t);
        tr=tr-median(tr);

        tr=tr/get_threshold(tr,1);
        sp=find(Result{y(1,i),j}.spike(t));
        [CS_list CS]=find_CS(tr-movmedian(tr,150,2),Result{y(1,i),j}.spike(t),15,5);
        Result{y(1,i),j}.Complex=CS>0;
        CS=CS>0; tr_CS=tr; tr_CS(~CS)=NaN;
        plot(t,tr+g*20,'color',cmap(i,:))
        hold all
        plot(sp,tr(sp)+g*20,'r.')
        plot(t,tr_CS+g*20,'color',[1 0 0.7])

        g=g+1;
        % ax1=[ax1 nexttile([1 1])];
        % plot(Result{y(1,i),j}.blue)
        % linkaxes([ax1],'x')
    end
end

%%
AvailResult=1-cellfun(@isempty,Result);
cmap=turbo(size(y,2));
clear Adap_data
figure;
for i=1:size(y,2) % Cellc 
    gg=1;
    Av=find(AvailResult(y(1,i),:));
    for j=[y(2,i) setdiff(Av,y(2,i))] % Sessions
        bw=bwlabel(Result{y(1,i),j}.blue>0);
        if max(bw)>3 
            delay=20; g=1; b=1;
            while g
                blue_seg=find(bw==b);
                if sum(Result{y(1,i),j}.spike(blue_seg(1):blue_seg(end)+delay))>0
                    g=0;
                    Adap_data{i,1}=Result{y(1,i),j}.blue(blue_seg(1));
                end
                b=b+1;
            end

            for b=max(bw)-5:max(bw)-1
                blue_seg=find(bw==b);
                Adap_data{i,2}(gg,1)=Result{y(1,i),j}.blue(blue_seg(1));
                Nss=sum(Result{y(1,i),j}.spike(blue_seg(1):blue_seg(1)+500));
                Ncs=sum(Result{y(1,i),j}.Complex(blue_seg(1):blue_seg(1)+500) & Result{y(1,i),j}.spike(blue_seg(1):blue_seg(1)+500));
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

clear
clc;
cd '/Volumes/BHL18TB_D1/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_OptopatchData_Arrangement.xlsx'], ...
    'Sheet1', 'B5:K104');

save_to='/Volumes/BHL18TB_D1/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
%%
nTau=[-50:70];
sub_time=[-11 -7; -22 -17; -11 -7; -11 -7;-11 -7];
Stim_time=[5 10; 4 8; 8 14; 8 12; 13 18];
DMDimg=[]; Subthimg=[]; STAtraces=[]; StimTAtraces=[]; STAmovies=[]; StimMovies=[]; Subthimg_Stim=[]; DMDcontour=[];
g=1;

for f=95:99
    load([fpath{f} '/Result.mat'])
    load(fullfile(fpath{f},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
    mov_res = SeeResiduals(mov_mc,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    [blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data);

    DMDimg(:,:,g)=blueDMDimg;
    DMDcontour{g}=bluePatt{1};
    nROI=size(Result.traces,1);

    tr_norm= Result.traces-movprc(Result.traces,100,30,2);
    tr_norm= tr_norm./get_threshold(tr_norm,1);
    NormTrace=Result.traces./get_threshold(Result.traces,1);
    Result.spike= find_spike_bh(tr_norm,5,1);
    bwBlue=bwlabel(Result.Blue);
    BlueOnset=find([0 Result.Blue(2:end)-Result.Blue(1:end-1)]==1);
    FirstSpike_Pulse=[];
    for j=1:max(bwBlue)
        Pulse_period=(bwBlue==j);
        FirstSpike_Pulse=[FirstSpike_Pulse find(Result.spike(1,:).*Pulse_period,1,'first')];
    end
    rmvSp=find(FirstSpike_Pulse+nTau(end)>size(Result.traces,2) | FirstSpike_Pulse+nTau(1)<1);
    FirstSpike_Pulse(rmvSp)=[];
    STAtraces(:,:,g)= squeeze(mean(reshape(Result.traces(:,FirstSpike_Pulse'+nTau),nROI,[],length(nTau)),2));
    StimTAtraces(:,:,g)= squeeze(mean(reshape(Result.traces(:,BlueOnset'+nTau),nROI,[],length(nTau)),2));

    TAmovie=squeeze(mean(reshape(mov_res(:,:,FirstSpike_Pulse'+nTau),sz(2),sz(1),[],length(nTau)),3));
    STAmovies(:,:,:,g)=-(TAmovie-mean(TAmovie(:,:,[1:7]),3));
    Subthimg(:,:,g) = mean(STAmovies(:,:,-nTau(1)+sub_time(g,1):-nTau(1)+sub_time(g,2),g),3);

    TAmovie=squeeze(mean(reshape(mov_res(:,:,BlueOnset'+nTau),sz(2),sz(1),[],length(nTau)),3));
    StimMovies(:,:,:,g)=-(TAmovie-mean(TAmovie(:,:,[1:7]),3));
    Subthimg_Stim(:,:,g) = mean(StimMovies(:,:,-nTau(1)+Stim_time(g,1):-nTau(1)+Stim_time(g,2),g),3);
    g=g+1;
end

SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.02,25);
interDendDist=[];
for i=1:nROI
    for j=1:nROI
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');
som_roi=find(dist_order==1);
geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)

%% Plot
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt);
Fimg=imgaussfilt(Result.ref_im,3);
Ref_filt=Result.ref_im;%-medfilt2(Result.ref_im,[20 20]);
ref_mask=Skeletonize_dendrite(Result.ref_im,10,0.02,25);
bound_im=zeros(size(ref_mask)); bound_im(bound:end-bound,bound:end-bound)=1;
ref_mask(~bound_im)=0; ref_mask=ref_mask>0;

figure(4); clf;
tiledlayout(5,2)
for n=1:5
    nexttile([1 1])
    show_im=imgaussfilt(Subthimg(:,:,n).*bound_im,3)./Fimg;%.*ref_mask;
    show_im=grs2rgb(show_im,turbo(100),0,0.02).*mat2gray(Ref_filt)*1.5;
    imshow2(show_im,[]); hold all
    plot(DMDcontour{n}([1:end 1],2),DMDcontour{n}([1:end 1],1),'color',[0 0.6 1],'LineWidth',1.5)
    title(['From ' num2str(-sub_time(n,1)) ' ms to ' num2str(-sub_time(n,2)) 'ms before spike'])
    colorbarHandle=colorbar; colormap(turbo);
    colorbarHandle.Ticks = [0 1];
    colorbarHandle.TickLabels = {'0', '2%'};
    title(colorbarHandle, '\DeltaF/F');
    nexttile([1 1])
    show_im=imgaussfilt(Subthimg_Stim(:,:,n).*bound_im,3)./Fimg;%.*ref_mask;
    show_im=grs2rgb(show_im,turbo(100),0,0.015).*mat2gray(Ref_filt)*1.5;
    imshow2(show_im,[]); hold all
    plot(DMDcontour{n}([1:end 1],2),DMDcontour{n}([1:end 1],1),'color',[0 0.6 1],'LineWidth',1.5)
    title(['From ' num2str(Stim_time(n,1)) ' ms to ' num2str(Stim_time(n,2)) 'ms after blue onset'])
    colorbarHandle=colorbar; colormap(turbo);
    colorbarHandle.Ticks = [0 1];
    colorbarHandle.TickLabels = {'0', '1.5%'};
    title(colorbarHandle, '\DeltaF/F');
end

figure(5); clf;
ROI_seq=[1 2 10 16 14];
offset=mean(STAtraces(:,1:15,:),2);
offset_stim=mean(StimTAtraces(:,1:15,:),2);
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt);
tiledlayout(5,3)
for n=1:5
    nexttile([1 1])
    l=plot(nTau,(STAtraces(dist_order,:,n)-offset(dist_order,:,n))'./F0(dist_order));
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2)); hold all
    ll=plot(nTau,squeeze(STAtraces(ROI_seq(n),:,n)-offset(ROI_seq(n),:,n))./F0(ROI_seq(n)),'Color','k','marker','*');
    title('Spike-triggered average')
    legend(ll,'Stim. ROI','Location','northwest')

    nexttile([1 1])
    imagesc(nTau,[1:nROI],(STAtraces(dist_order,:,n)-offset(dist_order,:,n))./F0(dist_order)',[0 0.04]); hold all
    plot(nTau(1)+1,find(ROI_seq(n)==dist_order),'r','Marker','>','MarkerFaceColor','r')
    colorbarHandle=colorbar; colormap(turbo);
    colorbarHandle.Ticks = [0 0.04];
    colorbarHandle.TickLabels = {'0', '4%'};
    title(colorbarHandle, '\DeltaF/F');

    nexttile([1 1])
    l=plot(nTau,(StimTAtraces(dist_order,:,n)-offset_stim(dist_order,:,n))'./F0(dist_order));
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2)); hold all
    ll=plot(nTau,squeeze(StimTAtraces(ROI_seq(n),:,n)-offset_stim(ROI_seq(n),:,n))./F0(ROI_seq(n)),'Color','k','marker','*');
    title('Stimulation-triggered average')
    legend(ll,'Stim. ROI','Location','northwest')

end

%% Low stim
for i=100
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg01.bin'],sz(2),sz(1)));
end
[blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data);
%% Low Stim plot
tr_norm= Result.traces-movprc(Result.traces,100,30,2);
tr_norm= tr_norm./get_threshold(tr_norm,1);
Result.spike= find_spike_bh(tr_norm,5,1);
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt);
NormTrace=Result.traces./F0';
subth_trace=get_subthreshold(NormTrace,Result.spike(1,:),5,10);

ROIinDMD=[4 5 6 8 12];

figure(6); clf;
tiledlayout(5,1)
nexttile([1 1])
imshow2(Result.ref_im,[]); hold all
plot(bluePatt{1}([1:end 1],2),bluePatt{1}([1:end 1],1),'color',[0 0.6 1],'LineWidth',2)
ax1=[];
title_str={'No Stim.','Low Stim.','No Stim.','Low Stim.'};
for n=1:4
ax1=[ax1 nexttile([1 1])];
t=[1:2000]+2000*(n-1);
imagesc(NormTrace(dist_order,t),[-2 4]*0.01); hold all
sp=find(Result.spike(1,t));
plot(sp,ones(length(sp),1)+nROI-1,'color','k','marker','^','MarkerFaceColor',[1 0 0],'LineStyle','none')
colormap(ax1(n),turbo);
title(title_str{n})
end

%% Grid Vertical Stim (soma)

for i=111
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
   % mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));
end
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(CamTrigger);
[blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data,'stack');
DMDtrigger=Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 3).data;
DMDtrigger=cumsum([0 DMDtrigger(2:end)-DMDtrigger(1:end-1)]==1)+1;
DMDtrigger=DMDtrigger(CamTrigger);
BlueOnset=[0 Blue(2:end)-Blue(1:end-1)]>0;
%%
figure(3); clf; cmap=turbo(size(bluePatt,2)); nTau=[-50:100];
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt);
tiledlayout(5,1)
nexttile([2 1])
show_footprnt_contour(blueDMDimg,Result.ref_im)
nexttile([3 1])
somtr=Result.traces(1,:);%-movmedian(Result.traces(2,:),200);
somtr=somtr./F0(1);
for d=1:max(DMDtrigger)
    t=find(DMDtrigger==d & BlueOnset);
    stimta=mean(somtr(t'+nTau),1);
    plot(nTau,stimta,'color',cmap(d,:)); hold all
end
axis tight

%% Grid Vertical Stim (dendrite)
for i=124
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
  %  mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));
end
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(CamTrigger);
[blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data,'stack');
DMDtrigger=Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 3).data;
DMDtrigger=cumsum([0 DMDtrigger(2:end)-DMDtrigger(1:end-1)]==1)+1;
DMDtrigger=DMDtrigger(CamTrigger);
BlueOnset=[0 Blue(2:end)-Blue(1:end-1)]>0;
%%
figure(3); clf; cmap=turbo(size(bluePatt,2)); nTau=[-50:100];
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt);
tiledlayout(5,1)
nexttile([2 1])
show_footprnt_contour(blueDMDimg,Result.ref_im)
nexttile([3 1])
somtr=Result.traces(5,:);%-movmedian(Result.traces(2,:),200);
somtr=somtr./F0(5);
for d=1:max(DMDtrigger)
    t=find(DMDtrigger==d & BlueOnset);
    stimta=mean(somtr(t'+nTau),1);
    plot(nTau,stimta,'color',cmap(d,:)); hold all
end
axis tight

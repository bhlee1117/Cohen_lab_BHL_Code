% Analysis on dendrite imaged from BU
% 2023/5/14, Byung Hun Lee

clear
[fpath] = uigetfile_n_dir;

%% Signal extraction

DAQ_rate=0.00001; isi_th=10;
for i=1:length(fpath)

disp([fpath{i}])
load([fpath{i} '/output_data.mat'])
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data

a=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
frm_rate=(a(2)-a(1))*DAQ_rate;
sz=Device_Data{1, 4}.ROI([2 4]);
mov_mc=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
mov_mc_g=arrayfun(@(z) imgaussfilt(mov_mc(:, :, z), 5), 1:size(mov_mc, 3), 'UniformOutput', false); %gaussian filter each frame
mov_mc_g = cat(3, mov_mc_g{:});

t=[1:size(mov_mc,3)-1]*frm_rate; t(t>length(Blue)*DAQ_rate)=[];
centers=[]; th=0.1;
[roi, ~]=clicky(mean(mov_mc,3));
mov_res=mov_mc_g-movmean(mov_mc_g,200,3);
centers=cell2mat(cellfun(@mean,roi,'UniformOutput',false)');

Result{i}.fp=fpath{i};
Result{i}.blue=Blue(round(t/DAQ_rate));
Result{i}.c_ftprnt=mask_footprint(centers,movmean(mov_res(:,:,1000:end),10,3),[],20);
Result{i}.traces=-(tovec(mov_mc_g)'*tovec(Result{i}.c_ftprnt))';
Result{i}.ref_im=mean(mov_mc,3);
Result{i}.spike=[];
tmp=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),100,2); tmp=tmp./get_threshold(tmp,1);
Result{i}.spike=find_spike_bh(tmp,4,3);
%tmp=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),100,2); tmp=tmp./get_threshold(tmp,1);
%Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3,1.5))>0;
for n=1:size(Result{i}.traces,1)
st=find(Result{i}.spike(n,:));
sp_seg=find(st(2:end)-st(1:end-1)>100);
sp_seg=[[1 sp_seg(1:end-1)+1]' sp_seg'];
s_locmax=zeros(1,size(Result{i}.spike,2));
for s=1:size(sp_seg,1)
    x=[sp_seg(s,1):sp_seg(s,2)];
    [~, arg]=max(Result{i}.traces(n,st(x)));
    s_locmax(st(x(arg)))=1;
end
Result{i}.STAmov{n}=generate_STAmov(s_locmax,mov_res,[10 40]);
end
%=xcorr(Result{i}.traces,tovec(mov_mc),5)


end


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

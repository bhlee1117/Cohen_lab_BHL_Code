% Analysis on AAV expression sample and plot, YQ601, ISO vs Ketamine
% 2023/08/27, Byung Hun Lee
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230901'
clear
[fpath] = uigetfile_n_dir;
%% Motion correction
for i=1:length(fpath)
    disp(['Motion correction processing on ' fpath{i}])
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 4}.ROI([2 4]));
    mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
    CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    mov_test=mov(:,:,250:350);

    try mov_test = single(mov_test)./single(max(mov_test.data(:)));
    catch disp('change to vm')
    mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end

    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));
    [mov_mc,xyField]=optical_flow_motion_correction_LBH(vm(mov(:,:,1:length(CamTrigger))),mov_ref, ...
        'normcorre');
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
    mcTrace=xyField;

    save([fpath{i} '/mcTrace.mat'],'mcTrace')
end

%% Signal extraction
for i=1:length(fpath)
    load([fpath{i} '/output_data.mat'])
    load([fpath{i} '/mcTrace.mat'],'mcTrace')
    sz=double(Device_Data{1, 4}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Encoder=Device_Data{1, 2}.buffered_tasks(1, 3).channels;
    CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    Blue=Blue(CamTrigger);
    nFrames=size(mov_mc,3);
    T_mean=squeeze(mean(mov_mc,[1 2]));
    avgImg=mean(mov_mc,3);
    bkTrace=get_blueoffTrace(T_mean,Blue,80);

    bkg = zeros(2, nFrames);
    bkg(1,:) = linspace(-1, 1, nFrames);  % linear term
    bkg(2,:) = linspace(-1, 1, nFrames).^2;  % quadratic term
    bkg(3,:) = bkTrace;

    mov_res= SeeResiduals(mov_mc,bkg,1);
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = -SeeResiduals(mov_res,mcTrace.^2);

    mov_mc_vec=tovec(mov_mc(:,:,1:nFrames)); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ref_im_vec=tovec(imgaussfilt(avgImg,1)); ref_im_vec=(ref_im_vec-mean(ref_im_vec))./std(ref_im_vec);

Result{i}.imcorr=sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1);    
c_ftprnt=mask_footprint(sz/2,mov_res,[],40);
Result{i}.c_ftprnt=imgaussfilt(c_ftprnt,2.5);
Result{i}.trace=[(tovec(mov_res)'*tovec(Result{i}.c_ftprnt))'];
Result{i}.Blue=Blue;
Result{i}.ref_im=avgImg;

end

save('Result_V2CheRiffST_20230902.mat','fpath','Result','-v7.3')

%%
clear cellnumber drug_ind
load('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230901/Result_V2CheRiffST_20230902.mat')
containString={'Awake','ket','CPP'};
for i=1:length(Result)
    sp=split(Result{i}.fpath,'_');
    cstr=find(contains(sp,'ell'));
    drug_str=find(contains(sp,containString));
    if isempty(drug_str)
        drug_ind(i)=0;
    else
        for j=1:length(containString)
        if ~isempty(find(contains(sp{drug_str},containString{j})))
        drug_ind(i)=j;
        end
        end
    end
    cellnumber(i)=str2num(cell2mat(regexp(sp{cstr},'\d+','match')));
end

Cell_list=unique(cellnumber);

for N=1:length(Cell_list)
Cell_ind=Cell_list(N);    
Scratch_list=find(cellnumber==Cell_ind);
[drug_sort, drug_sort_ind]=sort(drug_ind(Scratch_list),'ascend');
for d=1:4 % 4 states, Iso, Awake, Ket, CPP
Result_NList{N,d}=[];    
Scratch_state_list=find(drug_ind(Scratch_list)==(d-1));
for j=1:length(Scratch_state_list) %find cells match to state
n=Scratch_list(Scratch_state_list(j));
Result_NList{N,d}.trace(j,:)=Result{n}.trace;
tr=Result{n}.trace;
tr=tr-movmedian(tr,200); tr=tr./get_threshold(tr,1);
Result_NList{N,d}.Blue(j,:)=Result{n}.Blue; 
Result_NList{N,d}.spike(j,:)=find_spike_bh(tr,5,3);
end
end
end

%%

f1=figure; clf;
tiledlayout(5,1)
noi=11; Waveform=1;
DAQ_rate=0.00001;  
load([fpath{1} '/output_data.mat'])
legend_label={'Iso','Awake','Iso/Ket','Iso/CPP'};
CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));

cmap=distinguishable_colors(4);
t=CamTrigger*DAQ_rate;
g=1;
disp([num2str(Cell_list(noi)) 'th Cell'])
ax1=nexttile([4,1]);
for d=1:4
    if ~isempty(Result_NList{noi,d})
BlueShape=sum(Result_NList{noi,d}.Blue>0,2)>10000; %0 = pulse Ramp, 1=triangle
w=find(BlueShape==Waveform);
[~, arg]= sort(sum(Result_NList{noi,d}.spike(w,:),2),'descend');

    if ~isempty(w)
BlueWvf=rescale(Result_NList{noi,d}.Blue(w(arg(1)),:));        
plot(t,rescale(Result_NList{noi,d}.trace(w(arg(1)),:))+g,'color',cmap(d,:))
hold all
g=g+1;
    end
    end
end
ind=find(cellfun(@isempty,Result_NList(noi,:))==0);
legend(legend_label(ind),'location','northeastoutside')
axis off
ax2=nexttile([1,1]);    
plot(t,BlueWvf,'color',[0 0.6 1])
hold all
plot([t(end-1100) t(end-100)],[-1 -1],'color','k','linewidth',3)
text(mean([t(end-1100) t(end-100)]),-0.5,'1 sec','horizontalalignment','center')
axis off
linkaxes([ax1 ax2],'x')
set(f1, 'Position', [100, 100, 800, 400]);

saveas(f1,[num2str(noi) 'thCell_' num2str(Waveform) '.fig'])
print(f1, [num2str(noi) 'thCell_' num2str(Waveform) '.jpg'],'-djpeg', ['-r', num2str(600)]);




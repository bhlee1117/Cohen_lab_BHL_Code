% Analysis on AAV expression sample and plot, in house YQ201
% 2023/07/11, Byung Hun Lee

clear
[fpath] = uigetfile_n_dir;
%%
for i=9:length(fpath)
    load([fpath{i} '/output_data.mat'])
    load([fpath{i} '/mcTrace.mat'],'mcTrace')
    sz=double(Device_Data{1, 4}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Blue_shutter=Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 1).data;
    Blue=Blue.*Blue_shutter;
    Encoder=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    Encoder=movsum([0 abs(Encoder(2:end)-Encoder(1:end-1))],100);
    CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    Blue=Blue(CamTrigger);
    Encoder=Encoder(CamTrigger);

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
Result{i}.Encoder=Encoder;
Result{i}.ref_im=avgImg;

end

save('Result_V2CheRiffST_20230714.mat','fpath','Result','l','-v7.3')

%%
figure;
g=1;
for i=l
plot(rescale(Result{i}.trace)+g)
g=g+1;
hold all
end

%%
Result_onlyGood=Result;
Result_onlyGood(l)=[];
fpath_onlyGood=fpath;
fpath_onlyGood(l)=[];
%%
figure(3); clf;

for i=1:size(Result_onlyGood,2)
    sp=split(fpath_onlyGood{i},'_');
    cstr=find(contains(sp,'ell'));
    cellname(i)=str2num(cell2mat(regexp(sp{cstr},'\d+','match')));
end
[c_list c_arg]=unique(cellname);
cmap=distinguishable_colors(length(c_list));

for i=1:length(cellname)
plot(rescale(Result_onlyGood{i}.trace)+i,'color',cmap(find(c_list==cellname(i)),:))
hold all
plot(Result_onlyGood{i}.Blue/3+i,'color',[0 0.2 1])
end



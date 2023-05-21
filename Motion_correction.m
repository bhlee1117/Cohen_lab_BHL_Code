% Motion correction & Movie saving code



%% Load parameters
clear
[fname,fpath] = uigetfile('*.*','MultiSelect','on'); if fname == 0, return;end
frm_end=336000; f_seg=[[1:10000:frm_end] frm_end+1];

mov_test=vm([fpath],[1151:1250]);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3)); dim=size(mov_ref);
load([fpath '/settings.mat'])
%% Correction & Saving

for i=1:length(f_seg)-1
    i
    try
        mov=vm([fpath],[f_seg(i):f_seg(i+1)+10]);
    catch
        mov=vm([fpath],[f_seg(i):f_seg(i+1)-1]);
    end
    [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref);
    mcTrace{i} = squeeze(mean(xyField,[1 2]));

    mcTrace_hi = mcTrace{i}-movmean(mcTrace{i},50,1);
    mov_res = SeeResiduals(mov_mc,mcTrace{i});
    mov_res = SeeResiduals(mov_res,mcTrace_hi);
    mov_res = SeeResiduals(mov_res,mcTrace{i}.^2);

    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath 'motion_corrected/mc' num2str(i,'%02d') '.bin'])
    mov_res=vm(mov_res);
    mov_res.transpose.savebin([fpath 'motion_corrected/mc_res' num2str(i,'%02d') '.bin'])
end

save(fullfile(fpath,'Correction_result.mat'),'mcTrace','dim')

%%

im_G=imgaussfilt(mean(mov_res,3),2);
[centers radii]=Cell_segment_circle_080222(im_G);
centers=cell_detection_manual(mean(mov_res,3),centers);


c_ftprnt=mask_footprint(centers,mov_res,1:8000,6);
c_ftprnt=imgaussfilt(c_ftprnt,1);
%colr = flip(max(colormap(jet(size(c_ftprnt,3))),0),1);
colr = jet(size(c_ftprnt,3));
figure;imshow2(squeeze(sum(c_ftprnt(:,:,:).*reshape(colr,1,1,[],3),3)),[]);
traces = -tovec(c_ftprnt)' * tovec(mov_res);

%%
figure;
tiledlayout(6,4)
ax1 = nexttile([4 4]);
Blue=DAQ_waves.amplitude(4,round([1:size(traces,2)]*1.27*1e-3/1e-5));
plot(traces')
ax3 = nexttile([2 4]);
plot(Blue)
%%
volt=[]; volt_bin=[]; XY_traces=[];
load('20220830_voltage_trace.mat','c_ftprnt')
load('Correction_result.mat')
for i=1:length(f_seg)-1
    i
    try
        mov_mc=readBinMov_times([fpath 'motion_corrected/mc' num2str(i,'%02d') '.bin'],size(mov_test,1),size(mov_test,2),[1:10000]);
        XY_traces=[XY_traces; mcTrace{i}(1:10000,:)];
    catch
        mov_mc=readBinMov([fpath 'motion_corrected/mc' num2str(i,'%02d') '.bin'],size(mov_test,1),size(mov_test,2));
        XY_traces=[XY_traces; mcTrace{i}(:,:)];
    end
    mov_vec=double(tovec(mov_mc));
    tmp=tovec(c_ftprnt)'*mov_vec;
    tmp_bin=tovec(c_ftprnt>0.1)'*mov_vec;
    volt=[volt -(tmp-median(tmp,2))];
    volt_bin=[volt_bin -(tmp_bin-median(tmp_bin,2))];
end
save("20220915_voltage_trace_mc.mat","volt",'c_ftprnt');
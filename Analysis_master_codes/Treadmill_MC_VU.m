clear
clc;

[fpath] = uigetfile_n_dir();
%%
for i=1%:length(fpath)

    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[8000:9100];

    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:10000:frm_end] frm_end+1];

    if length(ref_time)>1000
        mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+1000]));
    else
        mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    end

    options_rigid = NoRMCorreSetParms('d1',size(mov_test,1),'d2',size(mov_test,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
    tic; [mov_test,shifts1,template1,options_rigid] = normcorre(mov_test,options_rigid); toc
    mov_test=vm(mov_test);
    mov_test = single(mov_test)./single(max(mov_test.data(:)));
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));

    for j=24:length(f_seg)-1
        g=1;
        try
            mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+10]));
        catch % when the image ends
            g=2;
            mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)-1]));
        end
        mov=vm(mov);
        [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

%         options_rigid = NoRMCorreSetParms('d1',size(mov_mc,1),'d2',size(mov_mc,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
%         tic; [mov_mc,shifts1,template1,options_rigid] = normcorre(mov_mc,options_rigid); toc

        ave_im=mean(mov_mc,3);
        mov_mc=vm(mov_mc);
        mov_mc.transpose.savebin([fpath{i} '/mc' num2str(j,'%02d') '.bin'])

        %mcTrace = squeeze(mean(xyField,[1 2]));
        mcTrace=xyField;
        save([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

        %  clear mov_mc mov
    end
end



% Run this after running 2nd session of Analysis_20250331_FigureMakingSeeSaw
nBinMov=150;
%% Signal extraction

for f=[foi(14:end)]%length(fpath)
    f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    NormalizedTrace=(Result.normTraces)./Result.F0_PCA;
    Subthres=get_subthreshold(NormalizedTrace(Dist_order{f}(noi_dist{f}),:),max(Result.spike(1,:),[],1)>0,7,17);
    [V D]=get_eigvector(Subthreshold{f});
    Sub_trace=V'*Subthres;

    bound=5;
    BoundROI=zeros(sz(2),sz(1));
    BoundROI(bound:end-bound,bound:end-bound)=1;
    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));
    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound))/std(ref_im_vec,0,1);
    %frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    f_seg_real=[f_seg(1:end-1)' f_seg(2:end)'];
    f_seg_real(1:end-1,2)=f_seg_real(1:end-1,2)+overlap;
    f_seg_real(2:end,1)=f_seg_real(2:end,1)+overlap+1;

    take_window=repmat([1 time_segment],length(f_seg)-1,1);
    take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
    take_window(end)=mod(f_seg(end),time_segment);
    take_window(take_window==0)=time_segment;
    Blue_on_Seg=unique(ceil(find(Result.Blue)/time_segment));

    BlueOffFrame=Result.Blue==0;
    RunTime=Result.VR(end,:)>0.002;
    StimOnLap=unique(Result.VR(8,:).*(Result.Blue>0));
    TrainingEnvLap=unique(Result.VR(8,:).*(Result.VR(2,:)==1));
    StimOfflap=setdiff(unique(Result.VR(8,:)),[StimOnLap TrainingEnvLap]);
    FrameOfInterest= BlueOffFrame & RunTime & ismember(Result.VR(8,:),StimOfflap);

    perispike_time=unique(find(Result.spike(1,:))'+[-5:10]);
    periblue_time=unique(find(Result.Blue>0)'+[-500:500]);
    t_fit= (ind2vec(size(Result.traces_bvMask,2),periblue_time,1)==0) & (ind2vec(size(Result.traces_bvMask,2),perispike_time,1)==0);
    [bleaching_fit] = expfitDM_2(find(t_fit)',-mean(Result.traces_bvMask(:,t_fit))',[1:size(Result.traces_bvMask,2)]',[60000]);

    SpikeHeight=Result.SpikeHeight_fit;
    Bin_position=ceil(Result.VR(5,:)/115*nBinMov);
    Bin_position(setdiff([1:size(Result.VR,2)],find(FrameOfInterest)))=NaN;
    Mov_PTA=zeros(sz(2)*sz(1),nBinMov); Mov_residue_PTA=zeros(sz(2)*sz(1),nBinMov);

    for j=1:length(f_seg)-1
        pos_vec=zeros(nBinMov,diff(f_seg_real(j,:))+1);
        for p=1:nBinMov
            pos_vec(p,:)=Bin_position(f_seg_real(j,1):f_seg_real(j,2))==p;
        end
        Prob_vec=Bin_position(f_seg_real(j,1):f_seg_real(j,2));
        validInd=find(~isnan(Prob_vec));

        if ~isempty(validInd)
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        %mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);

        bkg = zeros(1, size(mov_mc,3));
        bkg(1,:) = bleaching_fit(f_seg_real(j,1):f_seg_real(j,2));  % linear term

        mov_res= mov_mc-median(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);
        mov_res= get_subthreshold(tovec(mov_res),Result.spike(f_seg_real(j,1):f_seg_real(j,2)),7,17);
        mov_res= toimg(mov_res,sz(2),sz(1));
        mov_res_residue= SeeResiduals(mov_res,Sub_trace(1:NSeesawComponent(f),f_seg_real(j,1):f_seg_real(j,2)),1);
        
        mov_res=tovec(mov_res);
        mov_res_residue=tovec(mov_res_residue);

        mov_res= mov_res./SpikeHeight(f_seg_real(j,1):f_seg_real(j,2));
        mov_res_residue= mov_res_residue./SpikeHeight(f_seg_real(j,1):f_seg_real(j,2));

        Mov_PTA=Mov_PTA+mov_res(:,validInd)*pos_vec(:,validInd)'./sum(FrameOfInterest);
        Mov_residue_PTA=Mov_residue_PTA+mov_res_residue(:,validInd)*pos_vec(:,validInd)'./sum(FrameOfInterest);

        disp([num2str(j) ' th segment is stacked']);
        else
            disp([num2str(j) ' has no valid index']);
        end
    end

    Mov_residue_PTA=mat2gray(Mov_residue_PTA)*255;
    Mov_PTA=mat2gray(Mov_PTA)*255;
    Mov_PTA=toimg(Mov_PTA,sz(2),sz(1));
    Mov_PTA=vm(Mov_PTA);
    Mov_PTA.transpose.savebin([fpath{f} '/Place_triggered_movie.bin'])

    Mov_residue_PTA=toimg(Mov_residue_PTA,sz(2),sz(1));
    Mov_residue_PTA=vm(Mov_residue_PTA);
    Mov_residue_PTA.transpose.savebin([fpath{f} '/Place_triggered_movie_res.bin'])
end
% Analysis on AAV expression sample and plot, iGluSnFR3 in vivo imaging

clear
[fpath] = uigetfile_n_dir;
%%
for i=1:length(fpath)
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 4}.ROI([2 4]));
    mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
    %mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[1:30000]));
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
    %mcTrace = squeeze(mean(xyField,[1 2]));

    % shifts_r = squeeze(cat(3,xyField(:).shifts));
    % shifts_nr = cat(ndims(xyField(1).shifts)+1,xyField(:).shifts);
    % shifts_nr = reshape(shifts_nr,[],ndims(mov_mc)-1,size(mov_mc,3));
    % %shifts_x = squeeze(shifts_nr(:,2,:))';
    % %shifts_y = squeeze(shifts_nr(:,1,:))';
    % mcTrace=squeeze(shifts_nr)';
    mcTrace=xyField;

    save([fpath{i} '/mcTrace.mat'],'mcTrace')

mov_mc=double(mov_mc);    
nFrames=size(mov_mc,3);
bkg = zeros(1, nFrames);
%bkg(1,:) = linspace(-1, 1, nFrames);  % linear term
%bkg(1,:) = linspace(-1, 1, nFrames).^2;  % quadratic term
bkg(1,:) = squeeze(mean(mov_mc,[1 2]));


avgImg=mean(mov_mc,3);
mov_res=mov_mc-median(mov_mc,3);
mov_res=  SeeResiduals(mov_res,bkg);
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);

ref_imfilt=mat2gray(imgaussfilt(mean(mov_mc,3)-imgaussfilt(mean(mov_mc,3),40),1.5));
%dFF=mov_res;
clear mov_filt
for z=1:nFrames
mov_filt(:,:,z)=imgaussfilt(mov_mc(:,:,z),2);
end
dFF=(mov_mc-movmedian(mov_filt,500,3))./movmedian(mov_filt,500,3);
dFF_filt=imgaussfilt3(dFF,[1 1 5]);

    mov_mc_vec=tovec(mov_mc(:,:,1:nFrames)); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ref_im_vec=tovec(imgaussfilt(avgImg,1)); ref_im_vec=(ref_im_vec-mean(ref_im_vec))./std(ref_im_vec);

imcorr=sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1);    

gen_dff_movie([fpath{i}],ref_imfilt,dFF_filt)
writeMov_wTrace([fpath{i} ,'_imcorr'],dFF_filt,[],100,0.1,[0 0.1],[],imcorr)
end

%%
i=1; edgebound=10;
load([fpath{i} '/output_data.mat'])
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
mov_mc=mov_mc-median(mov_mc);
ROI=[200 100 300 300];
mov_mc=mov_mc(ROI(1):ROI(1)+ROI(3),ROI(2):ROI(2)+ROI(4),:);
mov_mc=mov_mc(edgebound:end-edgebound,edgebound:end-edgebound,1:1000);
MovFilt=imgaussfilt3(mov_mc,[2 2 1]);
resize_factor=0.5;

MovFilt_resz=imresize(MovFilt,resize_factor);
d=size(MovFilt_resz);
MovVec=tovec(MovFilt_resz);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);   
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

eigImgs = toimg(V,d(1),d(2));
figure;
for ind=1:10; nexttile([1 1]); imshow2(eigImgs(:,:,ind),[]); title(num2str(ind)); end
figure;
nexttile([1 1]);
for ind=1:10
plot(rescale(MovVec'*V(:,ind))+ind)
hold all
end
%%
pca_ind=[1 2 3];
sz2=size(mov_mc);
footprint = mat2gray(toimg(mean(abs(V(:,pca_ind)).*D(pca_ind)',2), d(1), d(2)));
reconMov=mat2gray(MovFilt.*imresize(footprint,sz2([2 1])));
dF_recon=reconMov-mean(reconMov,3);
moviefixsc(dF_recon,[0 0.1])


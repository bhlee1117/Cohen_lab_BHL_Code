clear
clc;
cd '/Volumes/cohen_lab-1/Lab/Labmembers/Byung Hun Lee/Data/20240120/'
fpath{1} = '/Volumes/cohen_lab-1/Lab/Labmembers/Byung Hun Lee/Data/20240120/224306BHLm119_Optopatch';
fpath{2} = '/Volumes/cohen_lab-1/Lab/Labmembers/Byung Hun Lee/Data/20240120/223634BHLm119_Optopatch';

%% Motion correction

for i=1:2

load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[9000:10000]; overlap=200;
time_segment=11000;

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

if length(ref_time)>2000
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
else
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
end
mov_test=rollingShutter_correction(mov_test,Device_Data{1, 3}.exposuretime,'fusion');
mov_test=mov_test(:,:,2:end);
[mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
mov_test=vm(mov_test);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));

for j=1:length(f_seg)-1
    try
        mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+overlap]));
    catch % when the image ends
        mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
    end

    mov=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');
    mov=vm(mov(:,:,2:end));
    if j==1
        mov=mov(:,:,[1 1:size(mov,3)]);
    end

    [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

    ave_im=mean(mov_mc,3);
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

    %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
    mcTrace=xyField; % Normcorre
    save([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

    %  clear mov_mc mov
end

end

%%
for i=2
disp(fpath{i})
DAQ_rate=0.000005;
load([fpath{i} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[2000:3000];
mov_test=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));
avgImg=mean(mov_test,3);
figure(3); clf;
imshow2(avgImg,[])
g=1; ROIpoly=[];
while g
    h = drawpolygon('Color','r');
    if size(h.Position,1)==1 %no more ROI
        g=0;
    else
        ROIpoly=[ROIpoly; {h.Position}];
    hold all
    plot(h.Position(:,1),h.Position(:,2))
    end
end
close(figure(3));
Result{i}.ROIpoly=ROIpoly;

load(fullfile(fpath{i},['/mcTrace' num2str(1,'%02d') '.mat']));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);

mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
Result{i}.ref_im=mean(mov_mc,3);
[clickyROI, original_trace]=clicky(mov_mc);

mov_res= mov_mc-mean(mov_mc,3);
mcTrace.xymean=movmean(mcTrace.xymean,3,2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean);
mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,3));

[SeeRes_trace]=apply_clicky(clickyROI,mov_res);
figure(3); clf;
plot(rescale(original_trace)); hold all;
plot(rescale(SeeRes_trace)+1);

n_comp=5;
mov_filt=imgaussfilt3(mov_res,[3 3 0.1]);
movVec=tovec(mov_filt);
Npoly=size(Result{i}.ROIpoly,1);
ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);

for p=1:Npoly %each ROIs
    mask(:,:,p) = poly2mask(Result{i}.ROIpoly{p}(:,1), Result{i}.ROIpoly{p}(:,2), sz(2), sz(1));
    pixelList=find(tovec(squeeze(mask(:,:,p))));
    subMov = movVec(pixelList,:);
    covMat = subMov*subMov';  % PCA within each region
    [V, D] = eig(covMat);
    D = diag(D); 
    D = D(end:-1:1);
    V = V(:,end:-1:1);
    vSign = sign(max(V) - max(-V));  % make the largest value always positive
    V = V.*vSign;
    coeff = mat2gray(mean(abs(V(:,1:n_comp)).*D(1:n_comp)',2));
    ftprnt(pixelList,p)=coeff;
end

Result{i}.ftprnt=toimg(ftprnt,sz(2),sz(1));
figure(4); clf;
imshow2(squeeze(sum(toimg(Result{i}.ftprnt,sz(2),sz(1)).*reshape(jet(Npoly),1,1,[],3),3)),[]);
1
end
%% Signal extraction
time_segment=11000; overlap=200;

for i=1:2
Result{i}.traces=[];
Result{i}.traces_res=[];
Result{i}.mcTrace=[];
Result{i}.im_corr=[];
bound=5;
ref_im_vec=tovec(Result{i}.ref_im(bound:end-bound,bound:end-bound));

load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
mc=[]; mov_mc=[];

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result{i}.Blue=Result{i}.Blue(CamTrigger);

for j=1:length(f_seg)-1
    j
    mov=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov=mov(:,:,[take_window(j,1):take_window(j,2)]);
        mc=[mc; mcTrace.xymean([take_window(j,1):take_window(j,2)],:)];

        mov_mc(:,:,end+1:end+size(mov,3))=mov;
end
mov_mc=mov_mc(:,:,2:end);

    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:)); 
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(1, size(mov_mc,3));
%     bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
%     bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Result{i}.Blue],30),3000);
    % mov_res = SeeResiduals(mov_res,mc);
    % mov_res = SeeResiduals(mov_res,mc.^2);
    % mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result{i}.traces=[Result{i}.traces -(tovec(mov_res)'*tovec(Result{i}.ftprnt))'];
    Result{i}.mcTrace=[Result{i}.mcTrace; mc];
    Result{i}.im_corr=[Result{i}.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation


Result{i}.traces=Result{i}.traces(:,1:length(CamTrigger)-1);
Result{i}.mcTrace=Result{i}.mcTrace(1:length(CamTrigger)-1,:);

end
save('BriefResult_SomaStim_20240121.mat','Result','fpath','-v7.3')
%% Clean up and norm

 for i=1:2
    Result{i}.normTrace=Result{i}.traces./get_threshold(Result{i}.traces,1);
    Result{i}.spike=find_spike_bh(Result{i}.normTrace-movmedian(Result{i}.normTrace,300,2),5,3);
 end
% 
 save('BriefResult_SomaStim_20240121.mat','Result','fpath','-v7.3')
%%
 load('BriefResult_SomaStim_20240121.mat')
%%
tiledlayout(6,2)
for i=1:2
    nexttile([2 2])
    l=plot(rescale2(Result{i}.normTrace,2)'+[1:size(Result{i}.normTrace,1)]);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(size(Result{i}.normTrace,1)),2))

    nexttile([1 2])
    plot(Result{i}.Blue,'color',[0 0.6 1])

    %nexttile([1 1])
    %imshow2(imfuse(mat2gray(max(Result{i}.ftprnt,[],3)),Result{i}.ref_im),[])
end


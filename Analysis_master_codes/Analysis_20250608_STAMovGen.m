if exist('bAPPropsMat{f,1}','var') & exist('fpath','var') & exist('PCresult','var')
    disp('bAPPropsMat{f,1} checked')
else
    error('run Analysis_20250428_FigureMakingbAPAmp.m first')
end

%% Load parameter and Result
[nROI nTime]=size(PCresult{f}.Subthreshold);
time_segment=15000; bound=5; overlap=200;
nTauPeak=[50 50];
Spikeorder2show=1;
load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result')
load(fullfile(fpath{f},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
frm_end=EndFrame(f);

%% Movie saving
disp([num2str(f) 'th file is processing'])

f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
f_seg_real=[f_seg(1:end-1)' f_seg(2:end)'];
f_seg_real(1:end-1,2)=f_seg_real(1:end-1,2)+overlap;
f_seg_real(2:end,1)=f_seg_real(2:end,1)+overlap+1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
take_window(take_window==0)=time_segment;

perispike_time=unique(find(PCresult{f}.SpikeMat(1,:))'+[-5:30]); perispike_time(perispike_time<1 | perispike_time>nTime)=[];
periblue_time=unique(find(PCresult{f}.BlueStim>0)'+[-5:30]); periblue_time(periblue_time<1 | periblue_time>nTime)=[];
t_fit= (ind2vec(nTime,periblue_time,1)==0) & (ind2vec(nTime,perispike_time,1)==0);
[bleaching_fit] = expfitDM_2(find(t_fit)',-mean(Result.traces_bvMask(:,t_fit))',[1:size(Result.traces_bvMask,2)]',[100000 10000]);
SpikeHeight=Result.SpikeHeight_fit;

fields = {'frameStacked','Transition'};
for i = 1:length(fields)
    STAinfo.(fields{i}) = []; % Clear the specified field for the given stype
end
STAinfo.StackedMovieN=0;

Spike2cat=find(bAPPropsMat{f,1}.IsNA==0 & bAPPropsMat{f,1}.IsBlue==0 & emptycell==0);
SpikeTimeVec=ind2vec(nTime,bAPPropsMat{f,1}.SpikeFrame(Spike2cat),1);
disp(['N = ' num2str(length(Spike2cat)) ' spikes are averaging'])
DirtTrace=sum(Result.dirtTrace>0,1)>0;
BlueTrace=Result.Blue>0;

Mov_PeakTA=[];
F0img=[]; c=1;
STAinfo.frameStacked{c}=[];

for j=1:length(f_seg)-1
    fprintf('Processing %2.0f/%2.0f movie chunk \n',j,length(f_seg)-1);
    [fInd]=find(SpikeTimeVec(:,[f_seg_real(j,1):f_seg_real(j,2)]));
    fIndVec=SpikeTimeVec(:,[f_seg_real(j,1):f_seg_real(j,2)]);
    %omitVec=max([fIndVec; Result.CStrace([f_seg_real(j):f_seg_real(j+1)]); DirtTrace([f_seg_real(j):f_seg_real(j+1)]); BlueTrace([f_seg_real(j):f_seg_real(j+1)])]);
    if ~isempty(fInd)
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        mc=movmean(mcTrace.xymean([take_window(j,1):take_window(j,2)],:),5,1);
        bkg = zeros(1, size(mov_mc,3));
        bkg(1,:) = bleaching_fit(f_seg_real(j,1):f_seg_real(j,2));  % bleaching regress out
        bkg=bkg./mean(bkg);

        mov_res= mov_mc-median(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res = SeeResiduals(mov_res,bkg,1);
        mov_res=tovec(mov_res);
        mov_res= mov_res./SpikeHeight(f_seg_real(j,1):f_seg_real(j,2));
        % if isempty(F0img)
        %     F0img=get_F0img(toimg(mov_res,sz(2),sz(1)));
        %     STAinfo.F0img=F0img;
        % end
        %mov_res=mov_res./tovec(F0img);
        %mov_res_baseline=get_subthreshold(mov_res,omitVec,15,200);

        [~, AddMov, AddSpikeTime]=get_STA(mov_res,fIndVec,nTauPeak(1),nTauPeak(2));
        AddSpikeTime=(AddSpikeTime)+f_seg_real(j)-1;
        Mov_PeakTA=cat(3,Mov_PeakTA,permute(AddMov,[1 3 2]));
        STAinfo.StackedMovieN=STAinfo.StackedMovieN+length(AddSpikeTime);
        STAinfo.frameStacked{c}=[STAinfo.frameStacked{c} AddSpikeTime];
        %STAinfo.SpikeType{Spikeorder2show}=[STAinfo.SpikeType{Spikeorder2show}; SpType(find(ismember(bAPPropsMat{f,1}(:,6),AddSpikeTime)),1)];

        MovtoWrite=vm(double(Mov_PeakTA)*10^12+6000);
        Movinfo=whos('MovtoWrite');
        if Movinfo.bytes > 4.0*10^9 | j==(length(f_seg)-1)
            MovtoWrite.transpose.savebin(fullfile(fpath{f},['SpikeTrigger_movie_' num2str(c) '.bin']))
            disp('Move on to the next bin')
            STAinfo.Transition{Spikeorder2show}(c)=AddSpikeTime(end);
            c=c+1;
            STAinfo.frameStacked{c}=[];
            Mov_PeakTA=[];
        end
    else
        disp([num2str(j) ' has no valid index']);
    end
end
disp('Stacking finished')
save(fullfile(fpath{f},'STAinfo'),'STAinfo')
%% Load SS

SS_list=PCresult{f}.SpikeClassMat(1,:);
[~, ST_spikeClassMat, SS_list_frame]=get_STA(PCresult{f}.SpikeClassMat,SS_list,nTauPeak(1),-1);
Isnearby=squeeze(sum(ST_spikeClassMat,[1 3]))>0;
IsolatedSS_frame=SS_list_frame(Isnearby==0);

SS2search=cellfun(@(x) find(ismember(x,IsolatedSS_frame)),STAinfo.frameStacked,'UniformOutput',false);
AlignMov=[];
for c=1:length(STAinfo.frameStacked)
    fprintf('Reading %2.0fth file out of %2.0f stacked \n',c,length(STAinfo.frameStacked));
    fname=['SpikeTrigger_movie_' num2str(c) '.bin'];
    Movreadsub=double(readBinMov_times(fullfile(fpath{f},fname),sz(2)*sz(1),sum(nTauPeak)+1,SS2search{c}))-6000;
    AlignMov=cat(3,AlignMov,Movreadsub);
end

STAmovieIsolatedSS=reshape(mean(AlignMov,3,'omitnan'),sz(2),sz(1),[]);
save(fullfile(fpath{f},'STAmovieIsolatedSS.mat'),'STAmovieIsolatedSS','nTauPeak','IsolatedSS_frame','-v7.3');

%% Load STA movies and save as mp4
Spikeorder2show=1;
alignedMovFN = ['SpikeTrigger_movie_SpOrder' num2str(Spikeorder2show)];
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN '*.bin']));
load(fullfile(fpath{f},'STAinfo.mat'))

AlignMov=[];
nReadmov=min([length(alignmovlist) 5]);
for l=1:nReadmov
AlignMov=cat(3,AlignMov,readBinMov(fullfile(fpath{f},alignmovlist(l).name),sz(2)*sz(1),length(nTau)));
end

STAmovie=-toimg(mean(double(AlignMov),3),sz(2),sz(1))-6000;
STAmovie=imgaussfilt(STAmovie,2);
cax_sub=[prctile(STAmovie(:),1) prctile(STAmovie(:),99.9)];
colorSTA=grs2rgb(tovec(STAmovie),colormap(turbo),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,sz(2),sz(1),range(nTau)+1,[]);
colorSTA=permute(colorSTA,[1 2 4 3]);
if isfield(Result,'Structure')
    STAMov_sub_Struc=colorSTA.*mat2gray(Result.Structure)*3;
else
    STAMov_sub_Struc=colorSTA.*mat2gray(Result.ref_im-100)*2;
end

figure(161); clf;
writeMov4d(fullfile(fpath{f},['SpikeTAmovie_' num2str(Spikeorder2show)]), ...
    STAMov_sub_Struc(:,:,:,-nTau(1)+[-50:30]),[1:81],10,1,cax_sub)
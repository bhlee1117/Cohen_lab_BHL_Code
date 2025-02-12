function [STA_movie_align_total StackedSpike] = get_STAMovie(fpath,filename,spikeVec,nTau_front,nTau_Back)
bound=6;
%f_tmp='/Volumes/BHL18TB_D1/20240218/134705BHLm117_FOV2_VR2';

nTau=[-nTau_front:nTau_Back];
alignmovlist=dir(fullfile(fpath,[filename '*.tiff']));
for l=1:length(alignmovlist)
    delete(fullfile(fpath,alignmovlist(l).name));
end

StackedSpike=[];
spikeframes{1}=find(spikeVec);

time_segment=15000;
fileInd=ceil(spikeframes{1}/time_segment);
frameInd=mod(spikeframes{1},time_segment);

load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
g=0; alignmov_ind=1;

STA_movie_align_total=[];
for j=unique((fileInd))
    j
    mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath '/mcTrace' num2str(j,'%02d') '.mat']);

    if j==length(f_seg)-1
        mc=mcTrace.xymean;
    else
        mov_mc=mov_mc(:,:,1:time_segment);
        mc=mcTrace.xymean(1:time_segment,:);
    end

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(2, size(mov_mc,3));
    bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);

    STA_movie_align=[];
    for k=find(fileInd==j)
        if frameInd(k)+nTau(1)>1 && frameInd(k)+nTau(end)<size(mov_res,3)
            mov_seg=mov_res(bound:end-bound,bound:end-bound,frameInd(k)+nTau);
            %mov_seg=mat2gray(mov_seg);
            STA_movie_align = cat(3,STA_movie_align, mov_seg);
            g=g+1;
            StackedSpike(1,g)=1;
            StackedSpike(2,g)=time_segment*(j-1)+frameInd(k);
        else
            STA_movie_align = cat(3,STA_movie_align, zeros(sz(2)-bound*2+1,sz(1)-bound*2+1,length(nTau)));
            g=g+1;
            StackedSpike(1,g)=0;
            StackedSpike(2,g)=NaN;
        end
    end
    STAinfo = whos('STA_movie_align_total');
    if STAinfo.bytes>2.0*10^9
        disp('STA size is too big, not concatenated but saved.')
    else
        STA_movie_align_total=cat(3,STA_movie_align_total,STA_movie_align);
    end
    if ~isempty(STA_movie_align)

        % write_tif_stack(STA_movie_align,fullfile(f_tmp,[alignedMovFN '_' num2str(alignmov_ind) '.tiff']))
        % alignmovlist=dir(fullfile(f_tmp,[alignedMovFN '*.tiff']));
        write_tif_stack(STA_movie_align,fullfile(fpath,[filename '_' num2str(alignmov_ind) '.tiff']))
        alignmovlist=dir(fullfile(fpath,[filename '*.tiff']));
        if alignmovlist(alignmov_ind).bytes > 2.0*10^9
            disp('Move on to 2nd tiff')
            alignmov_ind = alignmov_ind + 1;
        end
    end
end
end

function rois = uiget_rois(rec_frames,varargin)
    cur_dir = pwd;
    [fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
    if contains(fname,'.mat')
        rois = load(fullfile(fpath,fname));
        rois=rois.rois;
        return
    end
    
    mov = vm(fpath);
    figure;moviesc(mov)
    
    if ~isempty(varargin)
        start_frame=varargin{1};
        frame_duration=varargin{2};
        seg_size=varargin{3};
    else
        start_frame=[];
        frame_duration=[];
        seg_size=[];
    end
    
    if is_seg_illum(mov,rec_frames)
        [mov_new,start_frame,frame_duration,seg_size] = convert_seg_illum_mov(mov,rec_frames);
        figure;moviesc(mov_new)
        rois =  bwboundaries(get_roi_mask(mov_new,start_frame,frame_duration,seg_size));
        
    else
        rois =  bwboundaries(get_roi_mask(mov,start_frame,frame_duration,seg_size));
    end
    
    rois = cellfun(@(x) flip(x,2),rois,'uniformoutput',false);
    cd(fpath)
    save('rois.mat','rois')
    cd(cur_dir)
end

function answer = is_seg_illum(mov,rec_frames)
    im=mov(1:rec_frames).mean;
    im = mat2gray(im);
    T=graythresh(im);
    im=imbinarize(im,T);
    im1 = ones(size(im));
    if (im(:)'*im1(:))/numel(im)<1/4
        answer=true;
    else
        answer=false;
    end

end

function [mov_new,start_frame,frame_duration,seg_size] = convert_seg_illum_mov(mov,rec_frames)
n_seq = floor(size(mov,3)/rec_frames)+round(mod(size(mov,3),rec_frames)/rec_frames);
mov_new = vm(zeros(size(mov,1),size(mov,2),rec_frames));
for i=1:n_seq
im=mat2gray(mov((i-1)*rec_frames+1:i*rec_frames).mean);
im = imgaussfilt(adapthisteq(im,'cliplimit',0.005),10);
T=graythresh(im);
im=imbinarize(im,T(1));

figure(10);imshow(im)

[bd,~]=bwboundaries(im,'noholes');
% [~,rm_idx]=rmoutliers(cellfun(@length,bd));
% bd(rm_idx)=[];
n_seg=length(bd); seg_size=ceil(ceil(sqrt(numel(mov(:,:,1).data)/(n_seg*n_seq)))/10)*10;

    for j=1:ceil(numel(mov(1).data)/seg_size^2/n_seq)
        [x_c,y_c] = centroid(polyshape(bd{j}(:,2),bd{j}(:,1)));
        row = ceil(y_c/seg_size);
        col = ceil(x_c/seg_size);
        mov_new((row-1)*seg_size+1:min(row*seg_size,size(mov,1)),(col-1)*seg_size+1:min(col*seg_size,size(mov,2)),:)=...
        mov((row-1)*seg_size+1:min(size(mov,1),row*seg_size),(col-1)*seg_size+1:min(size(mov,2),col*seg_size),(i-1)*rec_frames+1:i*rec_frames);
    end
end
start_frame=1;frame_duration = rec_frames;
end
m_q = [63 14];

cur_dir = pwd;
    
    multi_file_opt = 'yes';
    mov_hcal = vm([]);
    n = m_q(1)+1;
    m = m_q(1);
    q = m_q(2);
    while strcmp(multi_file_opt,'yes')
        pause(.5)
        [fname,fpath]=uigetfile(matlabroot,'Please select hadamard calibration sequence','*.bin');
        [~,~,ext] = fileparts(fname);
        
        mov_hcal = vm(cat(4,mov_hcal.data,vm(fpath).data-100));
        
        uifig = uifigure;
        multi_file_opt = uiconfirm(uifig,'Add another calibration file?',...
            ['current number =' num2str(size(mov_hcal,4))],'options',{'yes','no'},'closeFcn', @(h,e) close(uifig));
    end
    mov_hcal = vm(sum(mov_hcal.data,4));
    mov_hcal = vm(single(mean(reshape(mov_hcal.data,mov_hcal.rows,mov_hcal.cols,n,[]),4)));
    
    
    hadtraces = hadamard_bincode_nopermutation(m)'*2-1;
%     hadtraces = hadamard_bincode(m)'*2-1;
    
    ccorr = mov_hcal*hadtraces;
    clear mov_hcal;
    ccorr = (ccorr)./(imgaussfilt(max(abs(ccorr)),2));    
    m = whos('ccorr');
    [userview, systemview] = memory;
    if systemview.PhysicalMemory.Available < m.bytes*16
        error 'not enough memory for Had63c14 calibration'
    end
    ccorr = ccorr.imresize(2,'bilinear');
%     ccorr = ccorr.blur(2).imresize(2,'bilinear');
    ccorr = cat(3,ccorr,-1.*ccorr);
    [~,ind] = max(ccorr.data,[],3);
    clear ccorr
    rpattern = ones(size(ind));
    rpattern(ind>size(hadtraces,2)) = -1;
    ind(ind>size(hadtraces,2)) = ind(ind>size(hadtraces,2)) - size(hadtraces,2);
    sind = sparse(1:numel(ind),ind,rpattern);
    vmind = vm(full(sind),size(ind));
    clear ind sind mvals
    scorr = vmind.imresize(.5,'bilinear');
%     scorr = vmind.imresize(.5,'bilinear').blur(.5);
    clear vmind
%%    
%     figure;moviesc(scorr)
    figure;imshow(scorr.std,[])
    
    figure;moviesc(scorr)
%% reconstruct image from Hadamard sequence
%     pause(1)
    [fpath]=uigetdir;
    file_prefix = 'hadamard_p';
%%
    f_list = ls([fpath '*' file_prefix '*']);
%% 
had_all = [];
    for ii = 2:size(f_list,1)
    cur_f = strsplit(f_list(ii,:));
    mov_h = vm(fullfile(fpath,cur_f{1}));

    mov_h = mov_h.correct_blank_marker-100;
    ref = mean(mov_h);
    had = im_lowpass(mean(mov_h.*scorr.ffzpad(mov_h.imsz))*2);
    had_all = cat(3,had_all,had);
    end
%     clear mov_h
%     figure;imshow(ref,[],'initialmagnification','fit');title('Wide Field')
%     figure;im=imshow(had,[],'initialmagnification','fit');title('Hadamard')
    


function im_f = im_lowpass(im_bin)
[nr,nc] = size(im_bin);
fft_im = fftshift(fft2(im_bin));
band_width = round(max(size(im_bin))/4);
f_ker = zeros(size(im_bin));
f_ker(1:band_width,1:band_width)=1;
f_ker(end-band_width:end,1:band_width)=1; 
f_ker(1:band_width,end-band_width:end)=1;
f_ker(end-band_width:end,end-band_width:end)=1;
f_ker = fftshift(f_ker);
im_f = abs(ifft2(ifftshift(fft_im.*imgaussfilt(f_ker,10))));

end
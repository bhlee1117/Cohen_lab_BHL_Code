function mov_filt= medfilt2_mov(mov,sz)
mov_filt=zeros(size(mov));
for i=1:size(mov,3)
    mov_filt(:,:,i)=medfilt2nan(mov(:,:,i),sz);
end
end
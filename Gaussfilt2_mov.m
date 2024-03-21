function mov_filt= Gaussfilt2_mov(mov,sz)
mov_filt=zeros(size(mov));
for i=1:size(mov,3)
    mov_filt(:,:,i)=imgaussfilt(mov(:,:,i),sz);
end
end
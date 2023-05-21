function out = getc2img(mov)

if isa(mov,'vm')
    mov = mov.data;
end
dmov = mov-repmat(mean(mov,3),[1 1 size(mov,3)]);

out = sqrt(mean(dmov(:,:,2:end).*dmov(:,:,1:end-1),3));

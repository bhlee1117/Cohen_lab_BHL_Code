function binnedmov = bin_movie(mov,bin)

if bin ==1;
    binnedmov =mov;
else
    [Ysize, Xsize,nfrm] = size(mov);
    tmov2 = squeeze(mean(reshape(mov,[bin,Ysize/bin,Xsize,nfrm])));

    tmov3 = permute(tmov2,[2 1 3]);

    binnedmov = permute(squeeze(mean(reshape(tmov3,[bin,Xsize/bin,Ysize/bin,nfrm]))),[2 1 3]);
end
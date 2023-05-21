fn= 'ROI_clip.tiff';
[mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = CellsortPCA('ROI_clip.tiff');
%%
[PCuse] = CellsortChoosePCs(fn, mixedfilters);
%%
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, .2, 20);
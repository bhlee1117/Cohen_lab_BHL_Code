function intens = imnotchfilt(intens,offlist,sigmalist,amplist)
%imnotchfilt	Spatial notch filter an image or movie.
% m = imnotchfilt(m, offrc, sigmarc, nhar) cancels the spatial frequency
% content of an image or movie at selected spatial frequencies. offrc is
% the offset or location for the 1st peak in units of cycles per pixel.
% sigmarc (optional if nhar is omitted) determines the width of the
% gaussian peak used for cancellation, and nhar (optional) is the number of
% harmonics at which this peak is repeated.
% 
%   2016 Vicente Parot
%   Cohen Lab - Harvard University
%
%     if ~exist('notch_mask','var')
        [nr, ~, ~] = size(intens);
        tr = (-2\nr:nr/2-1)';
%         ur = tr/nr;
    %%
    %{
        im = log(abs(fftshift(fft2(had(:,:,1)))));
        imagesc(...
            uc,ur,...
            im)
        axis image
    %}
    %%
        if ~exist('offlist','var') || isempty(offlist)
            error('no freqs')
        end
        if ~exist('sigmalist','var') || isempty(sigmalist)
            sigmalist = 10;
        end
        if numel(sigmalist)~=numel(offlist)
            sigmalist = sigmalist(1).*ones(size(offlist));
        end
        if ~exist('amplist','var') || isempty(amplist)
            amplist = 1;
        end
        if numel(amplist)~=numel(offlist)
            amplist = amplist(1).*ones(size(offlist));
        end
            
    %%
        notch_mask = ones(nr,1);
        for it = 1:numel(offlist)
%         for ir = -nhar:nhar
%             for ic = -nhar:nhar
%                 if ~ir && ~ic
%                     continue
%                 end
                notch_mask = notch_mask .* (1-amplist(it).*exp(-(abs(tr)-offlist(it)).^2/sigmalist(it)));
%             end
%         end
        end
%     end
    intens = real(ifft(bsxfun(@times,fft(intens),fftshift(notch_mask))));
%     had = had2;
    %
%     imshow(fftshift(log(abs(fft2(had2()-0*mean(mean(had2)))))),[])
%     im = log(abs(fftshift(fft2(had2(:,:,1)))));
%     imagesc(...
%         uc,ur,...
%         im)
%     axis image
%     imshow(notch_mask)


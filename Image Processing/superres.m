function out = superres(mov, kernel, ds, dt);
% function out = superres(mov, kernel, ds, dt);
% Make a superresolution movie.
% mov is a single-event movie (i.e. assuming already averaged over many
% spikes)
% Kernel is the waveform to fit to each pixel
% ds is the downsample factor.  An integer > 1.
% dt = time interval between frames in original movie (ms).
% Adam Cohen Feb. 7, 2012

[ysize xsize nframes] = size(mov);
avgimg = mean(mov,3);

dmov = mov - repmat(avgimg, [1 1 nframes]);  % work with the deviations from mean intensity
dmov = mat2gray(dmov);

slength = 300; % number of high resolution timesteps to search for peak

kernel = kernel - min(kernel);

[~,taumax] = max(kernel);  % Find the central peak and exactly one period of the kernel
L = watershed(kernel);  % Find the local basins of attraction
[~, kernstart] = min(kernel(L == L(taumax-1)));  % Find the index of the minimum in the basin before the peak
kernstart = kernstart + min(find(L == L(taumax-1)));
[~, kernstop] = min(kernel(L == L(taumax + 1)));  % Find the index of the minimum in the basin after the peak
kernstop = kernstop + taumax;
kernel(1:kernstart-1) = 0;  % Set the regions outside the central period to zero.
kernel(kernstop+1:end) = 0;
figure(35)
plot(kernel)
title('Kernel')


tau = 0:nframes-1;
tau2 = (0:1/ds:nframes-1/ds); % make a high resolution time axis
kernhr = interp1(tau, kernel, tau2, 'linear');  % make a high resolution kernel
kernhr(end-ds + 1:end) = 0;  % correct for the NaN's outside the interpolation region
% plot(tau2, kernhr, tau, kernel)

corrmov = zeros(ysize, xsize, slength);
[jnk, tshift] = max(kernel); 
for tau1 = 1:slength;
    tmp = circshift(kernhr', tau1-slength/2);
    tmp = reshape(tmp, ds, nframes);
    kernsamp = mean(tmp, 1);
%     plot(tau, kernel, tau, kernsamp)
%     pause(.01)
    kernsamp = reshape(kernsamp, 1, 1, nframes);
    corrmov(:,:,tau1) = mean(dmov.*repmat(kernsamp, [ysize, xsize]), 3);
end;
   
figure(32)
plot(1:slength, squeeze(mean(mean(corrmov))), 'r-');
title('Correlation vs lag')
[rmax, rindx] = max(corrmov, [], 3);

supermov = zeros(ysize, xsize, slength);
gausskern1 = fspecial('gaussian', [slength 1], 2);
ampimg = max(dmov, [], 3) - min(dmov, [], 3);
ampimg = mat2gray(ampimg);


for r = 1:ysize;
    for c = 1:xsize;
        supermov(r,c,:) = ampimg(r,c)*circshift(gausskern1, -slength/2 +rindx(r,c));
    end
end;
supermov = mat2gray(supermov);

map1 = colormap('hot');
map1(:,1) = map1(:,1) - .0417;
% map2 = map1;
% map2(:,3) = map1(:,1);
% map2(:,1) = map1(:,3);

% supercolmov = .5*repmat(ampimg, [1 1 3 slength]);
% % roffset = median(supermov(:));
% roffset = 0;
% for j = 1:slength;
%     supercolmov(:,:,:,j) = supercolmov(:,:,:,j) + grs2rgb(.9*(supermov(:,:,j)-.1), map1, 0, 1);
% end;
% 
% figure(1)
% clf
% % while(1)
%     for j = 1:slength;
%         imshow(supercolmov(:,:,:,j), 'InitialMagnification', 'fit')
%         text(1, 5, ['t = ' num2str(1000*j*dt/ds) ' ms'], 'Color', [.99 .99 .99], 'FontSize', 16)
% %         title(num2str(j))
%         pause(.01)
% %         M(j) = getframe(gca);
%     end;
% % end;
out = supermov;


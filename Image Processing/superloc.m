%% Super Temporal Localization - required inputs are kernel (kernel1), dspikemov,
% spikemovie average (spikemovav), weight image (weightimg_opto), and excite (if you
% are exciting optically. I will clean this up. For now, see 
% optopatch_v4.m for complete code. DRH


fig = 0;

fig = fig+1;
figure(fig);
plot(t(1:length(kernel1)),kernel1);
[val peak] = max(kernel1);

tbeg = 1;
tend = length(kernel1);
nt = tend-tbeg+1;

kernel2=kernel1(tbeg:tend); %select only around the peak

fig = fig+1;
figure(fig);
plot(t(1:length(kernel2)),kernel2')

kernel2 = kernel2-mean(kernel2);    %mean is zero

kernel2 = kernel2/sqrt(kernel2'*kernel2);   %normalize


dkernel = (kernel2(3:end)-kernel2(1:end-2))/2;
dkernel = [0; dkernel; 0];


dkernel'*kernel2;

dspikemov2 = dspikemov(:,:,tbeg:tend)- repmat(mean(dspikemov(:,:,tbeg:tend),3),[1 1 nt]);




% smoothing
smoothmask = zeros(5,5);
width = 0.75;
masksize = 5;

[X, Y] = meshgrid(1:masksize, 1:masksize);

smoothmask = exp(-((X-mean(X(:))).^2+(Y-mean(Y(:))).^2)/(2*width^2));

fig = fig+1;
figure(fig)
imshow(smoothmask,[],'Initialmagnification','fit')

smoothdspike = zeros(ysize,xsize,nt);
for q=1:nt
smoothdspike(:,:,q) = conv2(dspikemov2(:,:,q),smoothmask,'same');
end

smoothdspike = dspikemov2;  %Comment out if you want to SMOOTH

spikemov2D = reshape(smoothdspike,xsize*ysize,nt);
amp = spikemov2D*kernel2;

dkernelNorm = (dkernel'*dkernel);

endCorrection = 0.5*(kernel2(nt)^2-kernel2(1)^2);

delay1D = (amp*endCorrection-spikemov2D*dkernel)./(amp*dkernelNorm-(spikemov2D*dkernel)*(endCorrection));

amp = amp./(ones(size(delay1D))-delay1D*endCorrection);
amp = reshape(amp,ysize,xsize);

delay1D=delay1D*dt;
delay2D=reshape(delay1D,ysize,xsize);

cmin3 = mean(delay2D(:))-0.001*std(delay2D(:));
cmax3 = mean(delay2D(:))+0.001*std(delay2D(:));

fig = fig+1;
figure(fig)
imshow(delay2D, [cmin3 cmax3],'InitialMagnification','Fit')
% imshow([mat2gray(max(dspikemov2,[],3)-min(dspikemov,[],3)),mat2gray(amp)],[],'InitialMagnification','Fit')

% Calculate the RMS deviation of each pixel from its neighbors to the
% north, south, east, and west

compound = zeros(ysize,xsize,5);
shifting = [0,0;0,1;1,0;-1,0;0,-1];

for i = 1:5;
    compound(:,:,i) = circshift(delay2D,shifting(i,:));
end
errimg = squeeze(std(shiftdim(compound,2)));
% 
% errimg = (((delay2D - circshift(delay2D, [0,1])).^2 + (delay2D - circshift(delay2D, [0,-1])).^2 + ...
%          (delay2D - circshift(delay2D, [1,0])).^2 + (delay2D - circshift(delay2D, [-1,0])).^2)/4).^0.5;
errimg = errimg(2:end-1,2:end-1);
% Looking at the log of errimg is easier
fig = fig+1;
figure(fig)
imshow(log(errimg), [],'InitialMagnification','Fit')


% Get some statistics on different ROIs (clickroi is in the image
% processing folder on the server)
[meanval, stddev, npix, pixvals]=clickroi(errimg,weightimg_opto(2:end-1,2:end-1));
meanval


%Make a Super Localization Movie

subframeT = 0.1; %ms

initialT = -1.5; %ms
finalT = 1.5; %ms

times = [initialT:subframeT:finalT];

fig=fig+1;
figure(fig)
hist(mat2gray(weightimg_opto));

threshweight = im2bw(mat2gray(weightimg_opto(2:end-1,2:end-1)),0.05);
threshweight = threshweight.*weightimg_opto(2:end-1,2:end-1);

fig=fig+1;
figure(fig)
imshow([weightimg_opto(2:end-1,2:end-1), threshweight],[],'initialmagnification','fit')

tsmoothing = sum(sum(errimg.*threshweight))/sum(sum(threshweight))
manualtSmooth = 'no';
spaceSmooth = 'no';

superlocmov = zeros(ysize,xsize,length(times));
for q = 1:length(times)
%     superlocmov(:,:,q) = exp(-(delay2D-times(q)*ones(ysize,xsize)).^2./(2*errimg).^2); 
 superlocmov(:,:,q) = exp(-(delay2D-times(q)*ones(ysize,xsize)).^2/(2*tsmoothing).^2); 
end


% fig=fig+1;
% figure(fig);
% while 1
%     for j = 1:length(times)
%         imshow(superlocmov(:,:,j),[cmin4 cmax4], 'InitialMagnification', 'fit')
%         title([num2str(times(j)) ' ms'], 'FontSize', 14, 'color', [0.01 0.01 0.01])
% %         saveas(gca, ['./SuperLoc Movies/' datnum '_SuperLocAPmovie_' num2str((j)) '.tif']);  
%         pause(.1)
%     end;
% end;


superlocmov1 = superlocmov.*mat2gray(repmat(mat2gray(weightimg_opto), [1 1 length(times)]));
superlocColormov = zeros(ysize,xsize,3,length(times));
colorimg3 = repmat(mat2gray(spikemovav), [1 1 3]);
colorimg3(:,:,3) = colorimg3(:,:,3)+0.4*excite;
% colorimg3(:,:,2) = colorimg3(:,:,2)+0.1*excite;
% colorimg3(:,:,1) = colorimg3(:,:,1)+0.01*excite;



cmin4 = 0;
cmax4 = 0.4;

for j = 1:length(times);
    superlocColormov(:,:,:,j) = 0.8*colorimg3 + .5*grs2rgb(superlocmov1(:,:,j), colormap(hot), cmin4, cmax4);
end;

fig = fig+1;
figure(fig)
while(1)
for j = 1:length(times);
    imshow(superlocColormov(:,:,:,j), 'InitialMagnification', 'Fit')
%     text(1,2,[num2str(times(j)) ' ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])
%     saveas(gca, ['./SuperLoc Movies/' datnum '_SuperLocAPmovie_' num2str((j)) '.tif']); 
    pause(.1)
end;
end;

save([datnum '_SuperLoc.mat'],'superlocColormov','superlocmov1', 'delay2D','times','excite','amp','tsmoothing','manualtSmooth','spaceSmooth','peak','val')


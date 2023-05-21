% creat the binary image sequence for the DMD for efficient image
% segmentation.  Assume image is 1024 wide by 512 tall for simple binary
% processing.  We can take an ROI of this movie for actual illumination if
% appropriate.

nreps = 1;      % number of times to replicate the signal
mfs = 8;        % the minimum feature size in pixels
nbits = floor(log2(mfs));   % the number of unneeded bits because of the minimum feature size
nrow = 512;
ncol = 1024;
nbitr = ceil(log2(nrow)) - nbits;
nbitc = ceil(log2(ncol)) - nbits;
nbit = nbitr + nbitc;

% create a column vector of pixel addresses and convert to binary along the
% third dimension
C = (0:nrow-1)';
Cb = zeros(nrow,1,nbitr-nbits);
for i = nbitr:-1:1;
    Cb(:,:,nbitr+1-i) = floor(C/2^(i+nbits-1));
    C = C - Cb(:,:,nbitr+1-i)*2^(i+nbits-1);
end

% create a row vector of pixel addresses and convert to binary along the
% third dimension
R = 0:ncol-1;
Rb = zeros(1,ncol,nbitc);
for i = nbitc:-1:1
    Rb(:,:,nbitc+1-i) = floor(R/2^(i+nbits-1));
    R = R - Rb(:,:,nbitc+1-i)*2^(i+nbits-1);
end

Cm = repmat(logical(Cb),[1,ncol,1]);
Rm = repmat(logical(Rb),[nrow,1,1]);
SstripeTemp = cat(3,Rm,Cm);
Sstripe(:,:,1:3:3*nbit-2) = SstripeTemp;
Sstripe(:,:,2:3:3*nbit-1) = ~SstripeTemp;
Sstripe(:,:,3:3:3*nbit) = false;

temp1 = false(nrow,ncol);
c = 1;
r = 1;
ScheckTemp = false(nrow,ncol,nbit);
for i = 1:nbit
    temp2 = temp1;
    if mod(i,2) == 1
        temp1 = Rm(:,:,r);
        r = r + 1;
    else
        temp1 = Cm(:,:,c);
        c = c + 1;
    end
    ScheckTemp(:,:,i) = xor(temp1,temp2);
end
Scheck(:,:,1:3:3*nbit-2) = ScheckTemp;
Scheck(:,:,2:3:3*nbit-1) = ~ScheckTemp;
Scheck(:,:,3:3:3*nbit) = false;

S = cat(3,Sstripe,Scheck);
S = repmat(S,[1,1,nreps]);

%% now try some random binary images
% at a first glance, it seems that the pixels are too close together, so
% likely every neuron will be stimulated every time.  Now introduce a
% scale factor to factor to increase the fluctuation size.

% De-magnified DMD pixel size is ~10/3.  Cell body is on the 10-25 um
% length scale.  Optimal "pixel" size is probably something like half the 
% average intermolecular spacing, which will depend on transfection
% efficiency.

% first, try by making a smaller image and expanding

sf = 2*mfs;         % scaling factor
ff = 0.5;       % illuminated fill factor
Conj = 1;       % 1 if to include conjugate images, 0 if not
nrand = 2*nbit*nreps;     % number of random images
Srand = false(nrow,ncol,nrand);
if Conj
    for i = 1:nrand
        I = imresize((rand(round(nrow/sf),round(ncol/sf)) < ff),[nrow,ncol]);
        Srand(:,:,3*i-2) = I;
        Srand(:,:,3*i-1) = ~I;
        Srand(:,:,3*i) = false;
    end
else
    for i = 1:2*nrand
        Srand(:,:,i) = imresize((rand(round(nrow/sf),round(ncol/sf)) < ff),[nrow,ncol]);
    end
end
% for i = 1:nrand
%         imshow(Srand(:,:,i)); 
%     if i <= nrand
%         title(['imopen, pic ',int2str(i),', fill fraction = ',num2str(sum(sum(Srand(:,:,i)))/(nrow*ncol))]);
%     else
%         title(['imresize, pic',int2str(i),', fill fraction = ',num2str(sum(sum(Srand(:,:,i)))/(nrow*ncol))]);
%     end
%     pause(.5)
% end

S = cat(3,S,Srand);
figure(1)
for i = 1:size(S,3)
    imagesc(S(:,:,i)); axis image off; colormap gray; caxis([0 1]);
    title(int2str(i))
    pause(.2)
    Mov(i) = getframe;
end

movie2avi(Mov, 'excitation_seq.avi', 'fps', 3, 'compression', 'none')

%%
% embed in larger matrix for DMD
nrowD = 1080;   % number of rows in the DMD
ncolD = 1920;   % # columns in the DMD

SD = false(nrowD,ncolD,size(S,3));
r_offset = round((nrowD-nrow)/2);
c_offset = round((ncolD-ncol)/2)-87;
SD(r_offset:r_offset+nrow-1,c_offset:c_offset+ncol-1,:) = S;

fid = fopen('neuron_segment_stim_mov.bin', 'w');
fwrite(fid, permute(SD,[2 1 3]), 'ubit1','b');     % DMD reads by rows, Matlab saves by columns, big endian ordering
fclose(fid);
% 
% figure(1)
% for i = 1:size(SD,3)
%     imagesc(SD(:,:,i)); axis image;
%     title([int2str(i), ' of ', int2str(size(SD,3))])
%     pause(.1)
% end

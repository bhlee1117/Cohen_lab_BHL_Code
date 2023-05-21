function mov = ReadTifStackGen(fname)

% Pre-allocating the memory and passing the image info into imread
% significantly reduces read time.  I think it helps imread locate the
% images within the tiff stack.  Has comparable speeds to ImageJ and works
% with TIFF files save from anywhere.
%
% INPUTS
% fname: the name (and path if necessary) of the tiff stack.
%
% OUTPUT
% mov: the 3D array, still in uint16.  A conversion to double will often 
% be desired after reading in the movie

% Kit Werley, 8/2013

ImInfo = imfinfo(fname);
nrow = ImInfo(1).Height;
ncol = ImInfo(1).Width;
NumImages = length(ImInfo);

mov = zeros(nrow,ncol,NumImages,'uint16'); 
for i=1:NumImages
   mov(:,:,i) = imread(fname,'Index',i,'Info',ImInfo);
   if mod(i,100) == 0
       [int2str(i) ' frames read']
   end
end
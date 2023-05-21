% note! pixels and pad must be powers of 2
cm=1;
inch=2.54;
planes = 1 ;% Number of planes that will be used in calculating the 3-D image
lambda = 780e-7; %The wavelength of light
f = 5.824;  % real units (cm)
pad =2*256; 
pixels = 256; %Number of Pixels in the SLM
iterations=100; %iterations for block projection algorithm
gaussianBeamRadius = 0.41*inch;    % real units
slmsize=2.0;
% end of parameters section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check array size
if (floor(pad/2)*2~=pad)
    display ('pad must be a power of 2!');
    return
end
if (floor(pixels/2)*2~=pixels)
    display ('pixels must be a power of 2!');
    return
end
if (floor(768/pixels)~=768/pixels)
    display('768 must be a multiple of pixels!');
    return
end
if (pixels>768)
    display('pixels must be less than or equal to 768!');
    return
end
if (floor(planes/2)==planes/2)
    display('planes should be odd!');
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate pixel size in each plane and initialize matrices
fieldofview=slmsize/pixels*(pad-pixels)+slmsize;
pix = lambda*f/fieldofview;   %Unit size in focal plane
dpix = slmsize/pixels;        %Unit size in the SLM plane
[x y]=meshgrid(1:pad,1:pad);  %x and y are indices into the matrix
xiu=(x-pad/2.0-1)*pix;        %xiu and etau are real units for the focal plane
etau=(y-pad/2.0-1)*pix;
xu=(x-pad/2.0-1)*dpix;        %xu and yu are real units at the DOE
yu=(y-pad/2.0-1)*dpix;
% initial conditions:
% Ao(x, y) = electric field amplitude profile in the SLM plane
% Io(x, y) = intensity pattern at the focal plane along the optical axis
Ao=sqrt(exp(-2*(xu.^2+yu.^2)/gaussianBeamRadius^2));
% add in the aperture of the SLM
% so that the SLM fits into pixels x pixels in the center of pad x pad
slmaperture=((x>(pad-pixels)/2.0) & (x<=(pad-pixels)/2.0+pixels) & (y>(pad-pixels)/2.0) & (y<=(pad-pixels)/2.0+pixels) );
Ao=Ao.*slmaperture;
Ao=Ao/(sqrt(sum(sum(abs(Ao).^2)))+eps);
power=sum(sum(abs(fft2(Ao)).^2));
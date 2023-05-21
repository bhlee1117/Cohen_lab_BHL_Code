function im_filt = cutoff_filter(im0,type,freq)

type = lower(type);


% I=rgb2gray(RGB); % convert the image to grey 
A = fft2(im0); % compute FFT of the grey image
A1=fftshift(A); % frequency scaling
% Gaussian Filter Response Calculation
[M N]=size(A); % image size
R=freq; % filter size parameter 
X=0:N-1;
Y=0:M-1;
[X Y]=meshgrid(X,Y);
Cx=0.5*N;
Cy=0.5*M;

% Filtered image=ifft(filter response*fft(original image))
switch type
    case 'highpass'
        Lo=exp(-((X-Cx).^2+(Y-Cy).^2)./(2*R).^2);
        Hi=1-Lo; % High pass filter=1-low pass filter
        K=A1.*Hi;
        K1=ifftshift(K);
        im_filt=ifft2(K1);
        
        
    case 'lowpass'
        Lo=exp(-((X-Cx).^2+(Y-Cy).^2)./(2*R).^2);
        J=A1.*Lo;
        J1=ifftshift(J);
        im_filt=ifft2(J1);
        
    case 'bandpass'
        Bd= exp(-((X-Cx).^2+(Y-Cy).^2)./(2*R(2)).^2)-...
            exp(-((X-Cx).^2+(Y-Cy).^2)./(2*R(1)).^2);
        L = A1.*Bd;
        L1=ifftshift(L);
        im_filt=ifft2(L1);
end

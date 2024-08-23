function filtered_image=bandpass_2d(image,low_freq,high_freq)
% Read an image
gray_image = image;         % Convert to grayscale if it's a color image

% Get the image size
[M, N] = size(gray_image);

% Perform FFT
F = fft2(double(gray_image));

% Shift zero-frequency component to the center
Fshift = fftshift(F);

% Create a frequency grid
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2);

% Create the bandpass filter mask
bandpass_mask = (D >= low_freq) & (D <= high_freq);

% Apply the mask to the shifted FFT
Fshift_filtered = Fshift .* bandpass_mask;

% Shift the zero-frequency component back to the original position
F_filtered = ifftshift(Fshift_filtered);

% Perform inverse FFT to get the filtered image
filtered_image = ifft2(F_filtered);

% Take the real part of the image (remove imaginary component due to numerical errors)
filtered_image = real(filtered_image);

%Display the original and filtered images
% figure(3);
% subplot(1,2,1);
% imshow(gray_image, []);
% title('Original Image');
% 
% subplot(1,2,2);
% imshow(filtered_image, []);
% title('Bandpass Filtered Image');
end
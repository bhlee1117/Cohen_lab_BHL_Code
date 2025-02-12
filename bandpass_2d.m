function filtered_image = bandpass_2d(image, low_signal_radius, high_signal_radius)
% BANDPASS_2D applies a bandpass filter to a 3D image stack (frame-wise).
% 
% Inputs:
%   image             - a 3D array where each slice image(:,:,z) is a 2D image.
%   low_signal_radius - minimum feature radius (in pixels) to preserve.
%   high_signal_radius- maximum feature radius (in pixels) to preserve.
%
% This function removes features **smaller than low_signal_radius** and **larger than high_signal_radius**.
%
% Frequency estimation:
%   - A feature with radius R (pixels) has a dominant frequency of about **1/(2*R)** cycles per pixel.
%   - Cutoff frequencies are thus:
%       f_low  = 1/(2 * high_signal_radius)
%       f_high = 1/(2 * low_signal_radius)

for z = 1:size(image, 3)
    % Extract the z-th frame
    gray_image = image(:,:,z);
    
    % Get image dimensions
    [M, N] = size(gray_image);
    
    % Compute 2D FFT
    F = fft2(double(gray_image));
    Fshift = fftshift(F);  % Center zero-frequency component
    
    % Create frequency grid in cycles per pixel
    u = (-N/2 : N/2-1) / N;  % Horizontal frequency axis
    v = (-M/2 : M/2-1) / M;  % Vertical frequency axis
    [U, V] = meshgrid(u, v);
    D = sqrt(U.^2 + V.^2);  % Radial frequency in cycles per pixel
    
    % Compute frequency cutoffs based on radius in pixels
    f_low  = 1 / (2 * high_signal_radius);  % Large features cutoff
    f_high = 1 / (2 * low_signal_radius);   % Small features cutoff
    
    % Create the bandpass mask
    bandpass_mask = (D >= f_low) & (D <= f_high);
    
    % Apply the bandpass mask in Fourier space
    Fshift_filtered = Fshift .* bandpass_mask;
    
    % Inverse FFT to reconstruct the filtered image
    F_filtered = ifftshift(Fshift_filtered);
    filtered_image(:,:,z) = real(ifft2(F_filtered));
end

end

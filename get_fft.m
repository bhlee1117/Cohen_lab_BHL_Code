function [power_spectrum frequencies]=get_fft(trace,Fs)
%Fs = 1000;             % Sampling frequency (Hz)
t = 0:1/Fs:1-1/Fs;     % Time vector
%trace=mcTrace(:,1);
% Compute the FFT
N = length(trace);     % Number of samples
fft_result = fft(trace); % Compute the FFT

% Compute the power spectrum
power_spectrum = abs(fft_result).^2 / N;
power_spectrum=power_spectrum(1:N/2+1);
% Frequency axis for plotting
frequencies = Fs * (0:(N/2)) / N;

end
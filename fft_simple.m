function [power_spectrum, frequencies]=fft_simple(signal,Fs)


s=size(signal);
if s(1)>s(2)
    signal=signal';
end

% Calculate the power spectrum
for i=1:size(signal,1)
N = size(signal,2);                % Length of the signal
fr = (0:N-1) * (Fs/N);    % Frequency axis
fft_result = fft(signal(i,:));          % Compute the FFT
ps = abs(fft_result).^2 / N^2; % Calculate the power spectrum

frequencies(i,:)= fr(1:round(length(fr)/2));
power_spectrum(i,:)= ps(1:round(length(ps)/2));

end
end
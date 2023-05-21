function [P1 f]=fft_transient(X,Fs)
X=X-median(X);
%figure;
T = 1/Fs;             % Sampling period
L = length(X);             % Length of signal
t = (0:L-1)*T;
Y=fft(X);
%Y = fft(squeeze(sum(max(mov_bin_s2(:,:,100:3000),[],2),1)));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1)
xlabel('f (Hz)')
ylabel('Power')

end
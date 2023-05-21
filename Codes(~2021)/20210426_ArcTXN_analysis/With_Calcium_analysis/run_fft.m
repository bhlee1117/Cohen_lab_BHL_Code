function [f P1]=run_fft(time_vector,signal,sw)
if sw==0 %non uniform
L = size(time_vector,2);             % Length of signal
Y=nufft(signal,time_vector);

P2 = abs(Y/L);
P1 = P2(1:round(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = 1/(1/30)*(0:(L/2))/L;
else
L = size(time_vector,2);     
Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:round(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = 30*(0:(L/2))/L;
end
end
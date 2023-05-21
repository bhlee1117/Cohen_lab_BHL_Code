function plot_freq_spec(data_in,dt)

 L = length(data_in);
X = fft(data_in); 
P = abs(X/L); P = P(1:L/2+1,:); P(2:end-1,:) = 2*P(2:end-1,:);
F = 1/dt*(0:L/2)/L;

plot(F(2:end)',P(2:end,:))
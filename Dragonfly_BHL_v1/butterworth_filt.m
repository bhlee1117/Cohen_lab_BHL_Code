function data_filt = butterworth_filt(data,order,freq_range,rate)
if freq_range(1) ==0
    [b,a] = butter(order/2,freq_range(2)*2/rate,'low');
elseif freq_range(2) == inf
    [b,a] = butter(order/2,freq_range(1)*2/rate,'high');
else
    [b,a] = butter(order,freq_range*2/rate);
end
data_filt=filtfilt(b,a,data);
end
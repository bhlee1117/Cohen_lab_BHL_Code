function delay=get_delay(kymoTrace, binCounts)

% Function to detect the spike delay, spline interpolate and find the
% peak
% 2024.10.17 Byung Hun Lee, Cohen Lab
% Input: n (ROI) X T (peri-spike time),
%        binCounts, a number, number of interpolation bin: binCounts * T

x0=[1:size(kymoTrace,2)]-1;
xint=[x0(1):1/binCounts:x0(end)];
for s=1:size(kymoTrace,3)
for n=1:size(kymoTrace,1)

    if sum(isnan(kymoTrace(n,:,s)))>0
        delay(n,s)=NaN;
    else
kymoTrace_int = interp1(x0,kymoTrace(n,:,s),xint,'spline');
    [A,B] = max(kymoTrace_int);
    delay(n,s) = xint(B);
end
end
end
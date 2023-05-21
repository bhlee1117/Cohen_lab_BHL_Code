function out = subthreshfind(fluor, pulseperiod);
%subthreshold event peak finder
%DRH 03/26/12

cycles = length(fluor)/pulseperiod;
refluor = reshape(fluor, pulseperiod, cycles);
kernel = mean(refluor,2);
plot(circshift(kernel,-1))
[amp ind] = max(kernel);
nsearch = round(pulseperiod/10)*2;
amplitude = zeros(cycles,nsearch);
for i = 1:nsearch+1
    amplitude(:,i) = refluor'*circshift(kernel,-nsearch/2+i-1);
end

[amp1 ind1] = max(amplitude,[],2);
ind1 = ind1-nsearch/2+ind;
out = (0:pulseperiod:(cycles-1)*pulseperiod)+ind1';
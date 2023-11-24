function out = rem60hz(in, sampR, deltaF)
% function out = rem60hz(in, sampR, deltaF);
% Removes 60 Hz noise and harmonics from a data trace.
% in: input data trace
% sampR: sample rate of data trace
% deltaF: Frequencies within F +- deltaF/2 are set to zero.  Typically choose deltaF = 1.
% AEC 4 Sept. 2011

L = length(in);
freq = (0:(L-1))*sampR/L;
maxF = freq(round(L/2));

badfreqs = 60:60:maxF;
badindx = [];
for j = badfreqs;
    badindx = [badindx find((j-deltaF/2 < freq) & (j + deltaF/2 > freq))];
end;

infft = fft(in);
infft(badindx) = 0;
infft(L + 2 - badindx) = 0;

out = ifft(infft);

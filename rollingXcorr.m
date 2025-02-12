function [rollingcorr t_lag center_t]=rollingXcorr(VecX,VecY,windowSize,skipFrm)

if length(VecX)~=length(VecY)
    error('Check the length of Vectors')
end

nTime=length(VecX);
rollingcorr=[];
h = waitbar(0, 'Rolling Xcorr...');
g=1;
for t=floor(windowSize)+1:skipFrm:nTime-floor(windowSize)
t_region=[t-floor(windowSize):t+floor(windowSize)];
[rollingcorr(:,g) t_lag]=nanXCorr(VecX(1,t_region),VecY(1,t_region),windowSize, 0);

progress = t / (nTime-floor(windowSize/2)); % Calculate progress as a fraction
    waitbar(progress, h, sprintf('Progress: %d%%', round(progress * 100)));
    g=g+1;
end

center_t=floor(windowSize)+1:skipFrm:nTime-floor(windowSize);

close(h);
end

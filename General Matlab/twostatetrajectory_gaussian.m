function [odds,newmean0,newmean1,newvar0,newvar1,newk01,newk10] = twostatetrajectory_gaussian(data,mean0,mean1,var0,var1,k01,k10)
% funtion ile basladigin zaman matlab icin bir foksiyon olusturyorsun.
% oncelikel yapman 
odds = zeros(1,length(data));

stdratio = sqrt(var1/var0);

odds(1) = stdratio*exp((data(1)-mean1)^2/(2*var1) - (data(1)-mean0)^2/(2*var0));
for i = 2:length(data)
    odds(i) = stdratio*exp((data(i)-mean1)^2/(2*var1) - (data(i)-mean0)^2/(2*var0))*(odds(i-1)*(1-k01)+k10)/(odds(i-1)*k01+1-k10);
end

currtransitioncounts = zeros(2,2);
transitioncounts = zeros(2,2);

backwards = stdratio*exp((data(end)-mean1)^2/(2*var1) - (data(end)-mean0)^2/(2*var0));
statedatasum = data(end)/(1+odds(length(data)))*[odds(length(data));1];
statesqrdatasum = data(end)^2/(1+odds(length(data)))*[odds(length(data));1];

for i = length(data)-1:-1:1
    currtransitioncounts(1,1) = (1-k01)*odds(i)/(1+odds(i))*backwards/(1+backwards);
    currtransitioncounts(1,2) = k01*odds(i)/(1+odds(i))/(1+backwards);
    currtransitioncounts(2,1) = k10/(1+odds(i))*backwards/(1+backwards);
    currtransitioncounts(2,2) = (1-k10)/(1+odds(i))/(1+backwards);
    transitioncounts = transitioncounts + currtransitioncounts/sum(sum(currtransitioncounts));
    odds(i) = odds(i)*(backwards*(1-k01)+k01)/(backwards*k10+1-k10);
    statedatasum = statedatasum + data(i)/(1+odds(i))*[odds(i);1];
    statesqrdatasum = statesqrdatasum + data(i)^2/(1+odds(i))*[odds(i);1];
    backwards = stdratio*exp((data(i)-mean1)^2/(2*var1) - (data(i)-mean0)^2/(2*var0))*(backwards*(1-k01)+k01)/(backwards*k10+1-k10);
end

newmean0 = statedatasum(1)/(odds(length(data))/(1+odds(length(data)))+transitioncounts(1,1)+transitioncounts(1,2));
newmean1 = statedatasum(2)/(1/(1+odds(length(data)))+transitioncounts(2,1)+transitioncounts(2,2));
newvar0 = statesqrdatasum(1)/(odds(length(data))/(1+odds(length(data)))+transitioncounts(1,1)+transitioncounts(1,2)) - newmean0^2;
newvar1 = statesqrdatasum(2)/(1/(1+odds(length(data)))+transitioncounts(2,1)+transitioncounts(2,2)) - newmean1^2;
newk01 = transitioncounts(1,2)/(transitioncounts(1,1)+transitioncounts(1,2));
newk10 = transitioncounts(2,1)/(transitioncounts(2,1)+transitioncounts(2,2));

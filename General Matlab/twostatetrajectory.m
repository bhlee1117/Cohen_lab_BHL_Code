function [odds,newcountrate0,newcountrate1,newk01,newk10] = twostatetrajectory(data,countrate0,countrate1,k01,k10)
% funtion ile basladigin zaman matlab icin bir foksiyon olusturyorsun.
% oncelikel yapman 
odds = zeros(1,length(data));

photonodds = countrate0/countrate1;
expphotonratediff = exp((countrate1-countrate0));

k10dt = k10;
k01dt = k01;

odds(1) = expphotonratediff*photonodds^data(1);
for i = 2:length(data)
    odds(i) = expphotonratediff*photonodds^data(i)*(odds(i-1)*(1-k01dt)+k10dt)/(odds(i-1)*k01dt+1-k10dt);
end

currtransitioncounts = zeros(2,2);
transitioncounts = zeros(2,2);

backwards = expphotonratediff*photonodds^data(end);
statephotoncounts = data(end)/(1+odds(length(data)))*[odds(length(data));1];

for i = length(data)-1:-1:1
    currtransitioncounts(1,1) = (1-k01dt)*odds(i)/(1+odds(i))*backwards/(1+backwards);
    currtransitioncounts(1,2) = k01dt*odds(i)/(1+odds(i))/(1+backwards);
    currtransitioncounts(2,1) = k10dt/(1+odds(i))*backwards/(1+backwards);
    currtransitioncounts(2,2) = (1-k10dt)/(1+odds(i))/(1+backwards);
    transitioncounts = transitioncounts + currtransitioncounts/sum(sum(currtransitioncounts));
    odds(i) = odds(i)*(backwards*(1-k01dt)+k01dt)/(backwards*k10dt+1-k10dt);
    statephotoncounts = statephotoncounts + data(i)/(1+odds(i))*[odds(i);1];
    backwards = expphotonratediff*photonodds^data(i)*(backwards*(1-k01dt)+k01dt)/(backwards*k10dt+1-k10dt);
end

newcountrate0 = statephotoncounts(1)/((transitioncounts(1,1)+transitioncounts(1,2)));
newcountrate1 = statephotoncounts(2)/((transitioncounts(2,1)+transitioncounts(2,2)));
newk01 = transitioncounts(1,2)/(transitioncounts(1,1)+transitioncounts(1,2));
newk10 = transitioncounts(2,1)/(transitioncounts(2,1)+transitioncounts(2,2));

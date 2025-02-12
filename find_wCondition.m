function [location spikeInd]=find_wCondition(spike,cond)
% input: 
% spike = binary vector 1 x T
% cond = binary vector 1 x T
% find a spike within cond == 1 and give the index number of which spike is
% in condition
if length(spike) ~= length(cond)
    error('The length of spike and condition vector is different')
end
location=find((spike .* cond)>0);
bwspike=bwlabel(spike);
spikeInd=unique(bwspike.*cond);
spikeInd(spikeInd==0)=[];
end
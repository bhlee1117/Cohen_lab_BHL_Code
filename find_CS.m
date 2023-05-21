function [CS_list CS_spike]=find_CS(traces,spike,t_thresh,area_th)
if nargin<4
   area_th=0;
end
t_thresh=t_thresh/1.25; %ms
s_t=find(spike);
ISI=s_t(2:end)-s_t(1:end-1);
s_cluster=bwlabel(ISI<t_thresh);

CS_spike=zeros(1,size(traces,2));
CS_list=[];
for c=1:max(s_cluster)
s=find(s_cluster==c);

CS_list(c,:)=[s_t(s(1)) s_t(s(end)+1) mean(traces(s_t(s(1):s(end)+1)),'omitnan') length(s)+1]; %start, end, area, #of spike
end

if ~isempty(CS_list)
CS_list(find(CS_list(:,3)<area_th),:)=[];

for i=1:size(CS_list,1)
    CS_spike(1,CS_list(i,1):CS_list(i,2))=i;
end
end
end
%
function rm_ind=manual_trace_deletion(Result,noi)
seg=10000;
if nargin<2
noi=[1:size(Result.traces_res,1)];
end

rm_ind=[];
for i=1:ceil(size(Result.traces_res,2)/seg)
    g=1;
while g
r=Result.Arduino(:,4);
m=Result.mcTrace; m(~r,:)=NaN;
[f, ~]=show_traces_spikes(Result.traces_res(noi,:),Result.spike(noi,:),m);
xlim([seg*(i-1) seg*i])
%zoom('on','motion','horizontal')
zoom xon
[x y button]=ginput(1);
if isempty(button)
    g=0;
else
    hold all
plot([x x],[0 max(m(:))],'r','markersize',18)
[x2 y button]=ginput(1);
rm_ind=[rm_ind; round([x x2])];
end
%Result.traces_res(:,rm_ind(end,1):rm_ind(end,2))=NaN;
Result.spike(:,rm_ind(end,1):rm_ind(end,2))=0;
close(f)
end
end
end
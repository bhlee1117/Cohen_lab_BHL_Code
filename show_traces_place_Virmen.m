function show_traces_place_Virmen(traces,spike,place_bin,Virmen_data,poi,noi)
scale=-7;
vel_thresh=0;

Virmen_data(:,5)=Virmen_data(:,5)-min(Virmen_data(:,5))+1;
lap_length=max(Virmen_data(:,8));
Virmen_data(:,8)=round(Virmen_data(:,8));
%Ard_data(end-1:end,2:4)=repmat(Ard_data(end-2,2:4),2,1);
lap_end=[0; find(abs(Virmen_data(2:end,8)-Virmen_data(1:end-1,8))>0); size(Virmen_data,1)];
laps=[lap_end(1:end-1)+1 lap_end(2:end)];

for l=1:size(laps,1)
lap_trace(laps(l,1):laps(l,2))=l; end

cum_trace=Virmen_data(:,13);
vel_trace=(cum_trace(2:end)-cum_trace(1:end-1))';
vel_trace(end+1)=vel_trace(end);

lap_dist=max(Virmen_data(:,5));
bin_dist=ceil(Virmen_data(:,5)/((lap_dist)/place_bin));
cmap=jet(place_bin);
spike_run=spike; spike_run(:,vel_trace<vel_thresh)=NaN;
for p=1:place_bin
spot_stay=find(bin_dist==p); 
lap_spot_stay(p,:)=(bin_dist==p)'.*lap_trace;
end

traces_run=traces; traces_run(:,vel_trace<vel_thresh)=NaN;
S=traces; S(isnan(spike_run) | spike==0)=NaN;
%imagesc(lap_spot_stay)
figure;
ax1=[];
cmap=jet(max(lap_spot_stay(:)));
for p=poi
ax1=[ax1 nexttile([1 1])];
for l=1:max(lap_spot_stay(p,:))
    lines=plot(traces_run(noi,find(lap_spot_stay(p,:)==l))'+l*scale,'color',cmap(l,:));
    hold all
    plot(S(noi,find(lap_spot_stay(p,:)==l))'+l*scale,'r.','markersize',10);
end
title(['bin #' num2str(p)])
set(gca,'ytick',[1:length(noi)]*scale,'yticklabel',[noi])
end

Lap_FR=NaN(max(lap_spot_stay(:)),place_bin,length(noi));

for p=1:place_bin
    for l=1:max(lap_spot_stay(p,:))
        if ~isempty(find(lap_spot_stay(p,:)==l)) %only the laps the mouse went to the place bin
    Lap_FR(l,p,:)=mean(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(1*1e-3);
        end
    end
end

% nexttile([1 1])
% lines=plot([0 1],repmat([1 1],max(lap_spot_stay(:)),1)'+[1:max(lap_spot_stay(:))]);
% arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(flipud(cmap),2))
% text(zeros(1,size(laps,1)),[1:max(lap_spot_stay(:))]+1,num2str([max(lap_spot_stay(:)):-1:1]'))
% axis tight off
 linkaxes(ax1,'xy')

end

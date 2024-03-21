clear; clc;
cd '/Volumes/BHL_WD18TB/20231210_Arc_Calcium'
save_to='/Volumes/BHL_WD18TB/20231210_Arc_Calcium';
load(fullfile(save_to,"20210104_Data.mat"))

%% show traces

for f=1%:6

    f1=figure(10); clf;
f2=figure(11); clf;
    for d=1:3
velocity_threshold=2; place_bin=150;
valid_ind=find(cellfun(@(x) isempty(find(isnan(x))),dat{f}.Cal(:,d)));

Cal_trace=cell2mat(dat{f}.Cal(valid_ind,d));
Arc_onoff=dat{f}.Arc(valid_ind,d+1);
t_VR=Full_result{f}.VR{d}(:,1);
VR_trace=Full_result{f}.VR{d}(:,2)+2000; VR_trace(1)=0;
Lap_trace=zeros(1,length(VR_trace));
Lap_trace(find(abs(VR_trace(2:end)-VR_trace(1:end-1))>2000))=1;
Lap_trace=cumsum(Lap_trace)+1;
Cum_trace=VR_trace'+([1 Lap_trace(1:end-1)]-1)*3000;

isPC=cell2mat(dat{f}.isPC(valid_ind,d));
VR_data=[];
VR_data(1,:) = t_VR;
VR_data(5,:) = VR_trace;
VR_data(8,:) = Lap_trace;
VR_data(11,:) = Cum_trace;  
VR_data(12,:) = [0 Cum_trace(2:end)-Cum_trace(1:end-1)];  
sign={'off->on','on->off','on->on'};

for n=find(isPC)'

[Lap_FR Lap_V]=PlaceTrigger_average(Cal_trace(n,:),place_bin,VR_data,velocity_threshold,3000);
Lap_FR(find(sum(isnan(Lap_FR),2)>place_bin/2),:)=[];
if sum(ismember(Arc_onoff(n),[2 3 4]))>0
    f1=figure(10);
    
    nexttile([1 1])
imagesc(Lap_FR)
colormap(turbo)
    title([num2str(valid_ind(n)) 'th Cell is Arc ' char(sign(Arc_onoff(n)-1))])

else
    f2=figure(11);
    nexttile([1 1])
imagesc(Lap_FR)
colormap(turbo)
    title([num2str(valid_ind(n)) 'th Cell is Arc off'])
end
end
    end
    set(f1, 'Position', [100, 100, 1500, 600]);
saveas(f1,fullfile(save_to ,['Arc_on_file' num2str(f) '_day' num2str(d) '.fig']))
print(f1, fullfile(save_to ,['Arc_on_file' num2str(f) '_day' num2str(d) '.jpg']),'-djpeg', ['-r', num2str(400)]);

set(f2, 'Position', [100, 100, 1500, 600]);
saveas(f2,fullfile(save_to ,['Arc_off_file' num2str(f) '_day' num2str(d) '.fig']))
print(f2, fullfile(save_to ,['Arc_off_file' num2str(f) '_day' num2str(d) '.jpg']),'-djpeg', ['-r', num2str(400)]);

end

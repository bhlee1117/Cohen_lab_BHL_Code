function show_lap_rawtrace(volt,spike,trajectory,target_bin,noi,place_bin,lap_dist)

cmap=winter(place_bin);
bin_dist=ceil(trajectory/(lap_dist/place_bin));
%bin_dist(~run)=NaN;
l=bwlabel(ceil(trajectory/(lap_dist/place_bin))==place_bin)+1;
g=1; R=zeros(size(l,1),1);

for i=1:max(l); [a b]=max(find(l==i)); R(a)=1; end; R=cumsum(R);
for r=0:max(R)-1 %laps
    for p=1:place_bin %place bin
        tmp{r+1,p}=find(bin_dist==p & R==r);
    end
end

for n=noi
figure; rng=[4000 3000]; scale=get_threshold(volt(n,:),1)/range(volt(n,:))*10;
for i=1:max(R)-1
     
        plot([1:sum(rng)+1]*1.25,volt(n,tmp{i,target_bin}(1)-rng(1):tmp{i,target_bin}(1)+rng(2))+i*scale,'color',cmap(i+1,:))
        hold all
        line([1 sum(rng)+1]*1.25,[i*scale i*scale],'color',[0.6 0.6 0.6]);
        text([-400]*1.25,i*scale,[num2str((tmp{i,target_bin}(1)-rng(1))*1.25*1e-3,'%.1f') 's'])
        text([sum(rng)+50]*1.25,i*scale,[num2str((tmp{i,target_bin}(1)+rng(2))*1.25*1e-3,'%.1f') 's'])
        s_tmp=find(spike(n,tmp{i,target_bin}(1)-rng(1):tmp{i,target_bin}(1)+rng(2)));
        plot([s_tmp]*1.25,volt(n,tmp{i,target_bin}(1)-rng(1)-1+s_tmp)+i*scale,'r.')
    
end
yyaxis right
line([rng(1) rng(1)]*1.25,[0 1],'color','c')
axis off
end
end
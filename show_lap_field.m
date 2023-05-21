function [lap_place_volt tmp]= show_lap_field(volt,trajectory,run,noi,place_bin,lap_dist,reward_pos,im_on)
if nargin<8
    im_on=false;
end

if im_on
    figure;
    tiledlayout(ceil(length(noi)/4),4)
end

bin_dist=ceil(trajectory/(lap_dist/place_bin));
bin_dist(~run)=NaN;
l=bwlabel(ceil(trajectory/(lap_dist/place_bin))==place_bin)+1;
g=1; R=zeros(size(l,1),1);



for i=1:max(l); [a b]=max(find(l==i)); R(a)=1; end; R=cumsum(R);

for r=0:max(R)-1 %laps
    for p=1:place_bin %place bin
        tmp{r+1,p}=find(bin_dist==p & R==r);
    end
end


for n=noi
    lap_place_volt{n}=NaN(max(R)-1,place_bin);
    for r=1:max(R) %laps
        for p=1:place_bin %place bin
            if isempty(tmp{r,p})
            lap_place_volt{n}(r,p)=NaN;    
            else
            lap_place_volt{n}(r,p)=mean(volt(n,(tmp{r,p})),'omitnan');
            end
        end
    end
    % voltage trace
    if im_on
        ax1 = nexttile([1 1]);
        %figure;
        lap_place_volt{n}=reshape(movmean(tovec(lap_place_volt{n}'),3),place_bin,[])';
        lap_place_volt{n}=lap_place_volt{n}-min(lap_place_volt{n}(:));
        lap_place_volt{n}=lap_place_volt{n}./range(lap_place_volt{n}(:));
        s=mean(lap_place_volt{n},1,'omitnan'); s=s-min(s); s=s./range(s);
        lap_place_volt{n}=[lap_place_volt{n}; zeros(1,place_bin); s];
        imagesc(([lap_place_volt{n}]),[0.1 0.9]);
        hold on;
        line(repmat(ceil(reward_pos/(lap_dist/place_bin)),1,2),[0 max(R)+0.5],'color','r','linestyle',':');
        text(1,-1,[num2str(n) 'th neuron'],'Fontsize',12)
        axis tight off

    end
end
end
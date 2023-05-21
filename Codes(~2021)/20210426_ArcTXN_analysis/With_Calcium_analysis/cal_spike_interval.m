function I=cal_spike_interval(dat,groups,p_th)

for m=1:size(dat,2)
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m}{g,1}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
        end
        list{m}{g,1}=find(sum(list{m}{g,1},2)==size(groups{g},2)/2); % cells that match the conditions
    end
    for d1=1:3
        
        for g=1:size(groups,2)
        I{g}{m,d1}=[];    
        Peaks=dat{m}.Peak(list{m}{g,1}',d1); % Cells consist of Peak arrays
        
        for p=1:size(Peaks,1) 
            if size(Peaks{p,1},2)>1
            peak_th=Peaks{p,1}(find(Peaks{p,1}(:,1)>p_th),:);
        if size(peak_th,1)>1
            interval=peak_th(2:end,3)-peak_th(1:end-1,3);
        else
            interval=[];
        end
        I{g}{m,d1}=[I{g}{m,d1};interval];
            end
        end
        end
    end
end

end
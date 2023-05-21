function [Peaks Tr]=peak_finding(cal_sigma,C_df,ini_fin)
transient=zeros(1,size(cal_sigma,2));
g=1;
Tr=[];
ini_fin(find(ini_fin(:,1)>size(cal_sigma,2)),:)=[];
for i=1:size(ini_fin,1)
    if isempty(find(ini_fin(i,:)==0))
        if ini_fin(i,2)>size(cal_sigma,2)
            ini_fin(i,2)=size(cal_sigma,2);
        end
               
        [x arx]=max(cal_sigma(1,ini_fin(i,1):ini_fin(i,2)));
        arx=arx+ini_fin(i,1)-1;
        x=x/cal_sigma(1,ini_fin(i,2));
        t_min=0.2*log(x)/log(2);
        if 1%(ini_fin(i,2)-ini_fin(i,1)+1)*1/30>t_min
            transient(1,ini_fin(i,1):ini_fin(i,2))=1;
        Tr(g,1)=max(C_df(1,ini_fin(i,1):ini_fin(i,2)));  %max
        Tr(g,2)=size(C_df(1,ini_fin(i,1):ini_fin(i,2)),2); %length
        Tr(g,3)=sum(C_df(1,ini_fin(i,1):ini_fin(i,2))); % area
        g=g+1;
        else
            
        end
    else
    end
end
[pks,locs,widths,proms] = findpeaks(cal_sigma,[1:1:size(cal_sigma,2)]*1/30,'MinPeakProminence',2);
%[pks,locs,widths,proms] = findpeaks(cal_sigma,[1:1:size(cal_sigma,2)]*1/30);
locs=round(locs*30);
pks_cdf=C_df(1,locs);
if ~isempty(find(transient(locs)))
Peaks=[pks(find(transient(locs)))' pks_cdf(find(transient(locs)))' locs(find(transient(locs)))'*1/30 widths(find(transient(locs)))' proms(find(transient(locs)))' ];
else
Peaks=[];
end
end
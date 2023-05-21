function [dFF_area]=extract_area_dFF(dFF,ini_fin)
s=zeros(size(dFF,2),1);
if isempty(ini_fin)
dFF_area=0;
else
    if isnan(ini_fin)
        dFF_area=NaN;
    else
    g=1;
for i=1:size(ini_fin,1)
    if ini_fin(i,2)>size(dFF,2)
        if ini_fin(i,1)<size(dFF,2)
       s(ini_fin(i,1):end,1)=1; 
     
        end
    else
        %~isempty(dFF(1,ini_fin(i,1):ini_fin(i,2))) &&
        if  sum(ini_fin(i,:)==0)==0
        s(ini_fin(i,1):ini_fin(i,2),1)=1;
        g=g+1;
        end
    end
end
dFF(1,~s)=0;
dFF_area=sum(dFF);
end
end
end
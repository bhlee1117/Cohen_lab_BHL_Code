function [cal_sigma ini_fin]=cal_transient_detection(cal_trace,threshold,timescale)
plot([timescale:timescale:size(cal_trace,2)*timescale],cal_trace)
[x y]=ginput(2);
sd=std(full(cal_trace(round(x(1,1)):round(x(2,1)))));
m=mean(full(cal_trace(round(x(1,1)):round(x(2,1)))));

cal=(full(cal_trace)-m)/sd;
labelarray=bwlabel(cal>threshold);
ini_fin=[0 0];
g=1;
clear cal_tr
cal_tr=zeros(1,size(cal,2));
for i=1:max(labelarray)
    s=find(labelarray==i);
    if i>1
    if s(1,1)>ini_fin(g-1,2)
    ini_fin(g,1)=s(1,1);
    s=find(cal(1,ini_fin(g,1):end)<0.5);
    ini_fin(g,2)=ini_fin(g,1)+s(1,1)-1;
    cal_tr(1,ini_fin(g,1):ini_fin(g,2))=1;
    g=g+1;
    end
    else
    ini_fin(g,1)=s(1,1);
    s=find(cal(1,ini_fin(g,1):end)<0.5);
    ini_fin(g,2)=ini_fin(g,1)+s(1,1)-1;
    cal_tr(1,ini_fin(g,1):ini_fin(g,2))=1;
    g=g+1;    
    end
end
cal_sigma=[[timescale:timescale:timescale*size(cal,2)]; cal; cal_tr];
end
function [cal_sigma cal_transient ini_fin]=cal_transient_detection(ref_result,threshold,thresh_refrac,X,sig)
for i=1:size(X,1)
 
% m=mean(full(ref_result.C_or(i,round(X(i,1)):round(X(i,2)))));
% ref_result.C_df(i,:)=(ref_result.C_or(i,:)-m)/m;

sd=std(full(ref_result.C_df(i,round(X(i,1)):round(X(i,2)))));
m2=mean(full(ref_result.C_df(i,round(X(i,1)):round(X(i,2)))));
cal=sig*(full(ref_result.C_df(i,:))-m2)/sd; %chage a unit to sigma
labelarray=bwlabel(cal>threshold); % find evenets exceed threshold
ini_fin{i}=[0 0]; 
g=1;
cal_tr=zeros(1,size(cal,2));   
for j=1:max(labelarray) 
    s=find(labelarray==j);  % j th event
    if j>1 
    if s(1,1)>ini_fin{i}(g-1,2) % new event starts after previous evnet
    ini_fin{i}(g,1)=s(1,1);
    V=find(cal(1,ini_fin{i}(g,1):end)<thresh_refrac); %find when the event get back to baseline
    
if isempty(V) %the last transients
ini_fin{i}(g,2)=size(cal,2);
else
    ini_fin{i}(g,2)=ini_fin{i}(g,1)+V(1,1)-1;
end
    cal_tr(1,ini_fin{i}(g,1):ini_fin{i}(g,2))=1;
    g=g+1;
    end
    else % first 
   
%     try
    ini_fin{i}(g,1)=s(1,1);
    s=find(cal(1,ini_fin{i}(g,1):end)<thresh_refrac);
    if isempty(s)
    ini_fin{i}(g,2)=size(cal,2);    
    else
    ini_fin{i}(g,2)=ini_fin{i}(g,1)+s(1,1)-1;
    end
    cal_tr(1,ini_fin{i}(g,1):ini_fin{i}(g,2))=1;
    g=g+1;
%     catch
%     end
    end
end
cal_sigma(i,:)=cal;
cal_transient(i,:)=cal_tr;
clear cal_tr cal
end
end
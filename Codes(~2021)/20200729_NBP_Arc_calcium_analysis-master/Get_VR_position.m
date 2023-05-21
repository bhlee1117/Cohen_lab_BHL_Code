%function import_VR_pos
clear
frm_int=1/29.4809523809524;
save_name=['\\NEUROBIOPHYSICS\Virtual Reality\VR_image_data\JOE19'];
%%
format longG

for day=1:3
[fnm pth]=uigetfile('*.csv','Multiselect','on');
a=readmatrix([pth fnm]);
clear t
t = datetime(a(:,1),'ConvertFrom','datenum');
f=figure;
plot(t)
[t1 y]=ginput(1);
close(f);
f=figure;
plot(a(round(t1)-1500:round(t1)+1500,3))
[t2 y]=ginput(1);
close(f);
[m x]=min(a(round(t1)+round(t2)-1501-50:round(t1)+round(t2)-1501+50,3));
t_start=round(t1)+round(t2)-1501-51+x;
sw=1;
i=t_start;
clear rt
while sw 
rt(i,1)=(a(i,1)-(a(t_start,1)+0.00487))*(a(i+1,1)-(a(t_start,1)+0.00487)); 
if rt(i,1)<10^-11
    sw=0;
    t_end=i;
else
    i=i+1;
end
end
figure
try
plot(a(t_start-1000:t_end+1000,1),a(t_start-1000:t_end+1000,3))
catch
plot(a(t_start-1000:t_end,1),a(t_start-1000:t_end,3))
end
hold all
line([a(t_start,1) a(t_start,1)],[-3000 3000],'color',[1 0 0])
line([a(t_end,1) a(t_end,1)],[-3000 3000],'color',[1 0 0])
ylim([-3000 3000])
position{day}=[[frm_int:frm_int:frm_int*(t_end-t_start+1)]' a(t_start:t_end,3)];
end
save([save_name '.mat'],'position')
%%
function writeMov_plot(filename,ax1,nframe)
axis tight
%originalX=ax1.XLim;
%originalY=ax1.YLim;
originalX=[1 12000];
originalY=[0.5 36.5];
%filename='test.mp4';
%nframe=1000;
wind_size=[1500 37];
g=1;
list=[];
%set(ax1 ,'Layer', 'Top')
while g
    [x y]=ginput(1);
    if ~isempty(x)
        list=[list; repmat([x y],3,1)];
    else
        g=0;
    end
end

n=size(list,1)+1;
zoominLocation=[[list(:,1)-wind_size(1)/2 list(:,1)+wind_size(1)/2] [list(:,2)-wind_size(2)/2 list(:,2)+wind_size(2)/2]];
zoominLocation=[[originalX originalY]; zoominLocation ;[originalX originalY]];
lim_list=[[0:nframe/(n+2):nframe]' [[originalX originalY]; zoominLocation ;[originalX originalY]]];
lim_list_int=interp1(lim_list(:,1),lim_list(:,2:end),[1:nframe],'pchip');

myVideo = VideoWriter([filename],"MPEG-4"); %open video file 
myVideo.FrameRate = 70; %round(nframe/30);  %can adjust this, 5 - 10 works well for me
open(myVideo)


for i=1:nframe

%ax1.XLim=lim_list_int(i,[1 2]);
xlim(lim_list_int(i,[1 2]))
%ax1.YLim=lim_list_int(i,[3 4]);
ylim(originalY);
drawnow;
pause(0.005);
frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

close(myVideo)
ringBell()
end



    %end
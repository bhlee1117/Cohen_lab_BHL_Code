function place_cell_movie(filename,Arduino_data,volt,spike,noi,frame_down,mov_test,ref_mask)
cmap=autumn(length(noi));
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[0 0 0],'Renderer','painters','position',[100 100 1000 700]);

tread_im=imread('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Group Meeting Presentation/20220715_GroupMeeting_figures/Picture1.png');
tread_im=imrotate(tread_im,180);
total_frame=size(volt,2);
%frame_down=15000;
scale=5;
window=900;
t=[1:round(total_frame/frame_down):total_frame-window];

myVideo = VideoWriter(filename,"MPEG-4"); %open video file
myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
%v=volt-median(volt,2);
tmp=volt; tmp(spike==1)=NaN;
v=(volt-median(volt,2))./std(tmp,0,2,'omitnan');
open(myVideo)

for i=t
    clf;
    handles.axes1 = axes('Units','pixels','Position',[50 450 250 200],'Color',[0 0 0]);
    imshow2(mov_test.mean,[])
    hold on
    for n=1:length(noi)
        plot(ref_mask{noi(n)}(:,1),ref_mask{noi(n)}(:,2),'color',cmap(n,:),'linewidth',2)
        text(mean(ref_mask{noi(n)}(:,1))+25,mean(ref_mask{noi(n)}(:,2)),num2str(n),"HorizontalAlignment","left",'color',cmap(n,:),'FontSize',15)
    end
    axis off


    handles.axes1 = axes('Units','pixels','Position',[350 400 630 100],'Color',[0 0 0]);
    imagesc(tread_im)
    axis tight equal off
    hold on
    plot(Arduino_data(i,2)/max(Arduino_data(:,2))*size(tread_im,2),size(tread_im,1)/2,'r.','markersize',24)

    handles.axes2 = axes('Units','pixels','Position',[350 500 630 170],'Color',[0 0 0]);
    plot(Arduino_data([1:i],1),Arduino_data([1:i],3)*4000,'color',[0 0.6 1],'linewidth',1.5)
    hold all
    plot(Arduino_data([1:i],1),Arduino_data([1:i],2),'color','w','linewidth',2.5);
    
    
    ylim([0 max(Arduino_data(:,2))])
    xlim([0 max(Arduino_data(:,1))])
  axis off

    handles.axes3 = axes('Units','pixels','Position',[50 10 940 420],'Color',[0 0 0]);
    for n=1:length(noi)
        
            if i<=window
                plot([NaN(1,2*window+1-(i+window)) v(noi(n),1:i+window)]+scale*n,'color',cmap(n,:))
                hold all
                plot([[1:2*window+1-(i+window)] 2*window+1-(i+window)+find(spike(noi(n),1:i+window))],[NaN(1,2*window+1-(i+window)) v(noi(n),find(spike(noi(n),1:i+window)))]+scale*n,'wo')
                text(-20,scale*n,num2str(n),"HorizontalAlignment","left",'color',cmap(n,:),'FontSize',15)
            else
                plot(v(noi(n),i-window:i+window)+scale*n,'color',cmap(n,:))
                hold all
                plot(find(spike(noi(n),i-window:i+window)),v(noi(n),find(spike(noi(n),i-window:i+window))+i-window-1)+scale*n,'wo')
                text(-20,scale*n,num2str(n),"HorizontalAlignment","left",'color',cmap(n,:),'FontSize',15)
            end

 

    end
    line([window window],[0 scale*(length(noi)+1.7)],'color','w','linestyle','--')
    text([-10], [scale*(length(noi)+2)],[num2str(Arduino_data(i,1),3) ' sec'],'color',[1 1 1])
       xlim([0 window*2+1])
       ylim([0 scale*(length(noi)+1.7)])
    axis off
    %title([num2str(i/frmrate*1000) ' frame'])
    pause(0.005) %Pause and grab frame

    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)
end

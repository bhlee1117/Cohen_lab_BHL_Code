function gen_dff_movie(filename,ref_im,mov_im,Dwave)
if nargin<4
    Dwave=[];
end
%ref_im=mean(mov_mc,3); 
ref_im=ref_im./max(ref_im(:));
%mov_im=-movmean((mov_mc-movmean(mov_mc,150,3)),20,3);
mov_im=mov_im./reshape(range(tovec(mov_im),1),1,1,[]);
mov_im(mov_im<0)=0; 
mov_im=mov_im*1/prctile(mov_im(:),99.9);
cmap=[0.4 0.4 0.4; 1 0 0];
%%
mov_merge=zeros(size(ref_im,1),size(ref_im,2),size(mov_im,3),3);

for c=1:3
    mov_merge(:,:,:,c)=mov_merge(:,:,:,c)+repmat(ref_im,1,1,size(mov_im,3))*cmap(1,c);
    mov_merge(:,:,:,c)=mov_merge(:,:,:,c)+mov_im*cmap(2,c);
end

%%

f=figure;
myVideo = VideoWriter(filename,"MPEG-4"); %open video file
myVideo.FrameRate = 300; %round(nframe/30);  %can adjust this, 5 - 10 works well for me
open(myVideo)
for i=1:size(mov_merge,3)
    clf;
    imagesc(squeeze(mov_merge(:,:,i,:)),[0 3])
    axis equal tight off
    title([num2str(i/800,'%2.2f') ' sec'])
    hold all;
     if ~isempty(Dwave)
text(5,-12,['blue int:' num2str(Dwave(i))],'color',[0 0.6 1],'FontSize',12,'HorizontalAlignment','left')
if Dwave(i)
    plot(15,15,'.','color',[0 0.6 1],'markersize',32)
end
    end
    drawnow;
    pause(0.005);
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

close(myVideo)
close(f)
ringBell()
end
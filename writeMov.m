function writeMov(filename,movie,stack,rate,frmrate,range,dmd_mask_sequence_rois,Dwave)
if nargin < 7;
    dmd_mask_sequence_rois=[];
    Dwave=[];
end


if isempty(stack)
    stack=[1:size(movie,3)];
end

if isempty(range)
    range=[min(tovec(movie(:,:,floor(size(movie,3)/2)))) max(tovec(movie(:,:,floor(size(movie,3)/2))))]
end

figure(1);

myVideo = VideoWriter([filename],"MPEG-4"); %open video file 
%myVideo = VideoWriter([filename],"Uncompressed AVI");
myVideo.FrameRate = rate;  %can adjust this, 5 - 10 works well for me
myVideo.Quality= 100;
open(myVideo)

try
for i=stack
    clf;
    
    pbaspect([size(double(movie),2) size(double(movie),1) 1]*2)
    imshow2(movie(:,:,i),range)
    axis tight off equal
    hold on
    if ~isempty(dmd_mask_sequence_rois)
        %plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'color',[0 0.6 1],'linewidth',2)
        plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'ro','markersize',12)
    end
    if ~isempty(Dwave)
text(5,-12,['blue int:' num2str(Dwave(i))],'color',[0 0.6 1],'FontSize',12,'HorizontalAlignment','left')
if Dwave(i)
    plot(15,15,'.','color',[0 0.6 1],'markersize',32)
end
    end
    title([num2str(i/frmrate) ' ms'])
    %set(gca, 'Position', [100, 100, 1700, 800]);
    %title([num2str(i*frmrate) ' \mum'],'HorizontalAlignment','left')
    pause(0.005) %Pause and grab frame
%colormap(turbo)
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    %pause(0.1)
end
end
close(myVideo);
close(figure(1));
end
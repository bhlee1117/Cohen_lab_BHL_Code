function writeMov_wKymo(filename,movie,stack,rate,frmrate,range,dmd_mask_sequence_rois,Kymo)
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
open(myVideo)

Kymo=Kymo(:,stack);
try
for i=1:length(stack)
    clf;
    tiledlayout(4,1)
    ax1=nexttile([3 1]);
    imshow2(movie(:,:,stack(i)),range)
    colormap(ax1,'turbo')
    %colormap(ax1,Aurora)
    axis tight off equal
    hold on
    if ~isempty(dmd_mask_sequence_rois)
        plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'color',[0 0.6 1],'linewidth',2)
        %plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'ro','markersize',12)
    end
    if ~isempty(Kymo)
    ax2=nexttile([1 1]);
    imagesc(Kymo)
    hold all
    plot([i i],[0.5 size(Kymo,1)+0.5],'color','r')
    axis tight
    end
    colormap(ax2,turbo)
    title([num2str(stack(i)/frmrate) ' ms'])
    %title([num2str(i*frmrate) '\mum'],'HorizontalAlignment','left')
    pause(0.005) %Pause and grab frame

    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    %pause(0.1)
end
end
close(myVideo);
close(figure(1));
end
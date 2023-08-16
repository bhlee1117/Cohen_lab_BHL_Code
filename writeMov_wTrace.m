function writeMov_wTrace(filename,movie,stack,rate,frmrate,range,dmd_mask_sequence_rois,Trace)
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
try
for i=stack
    clf;
    tiledlayout(4,1)
    nexttile([3 1])
    imshow2(movie(:,:,i),range)
    axis tight off equal
    hold on
    if ~isempty(dmd_mask_sequence_rois)
        %plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'color',[0 0.6 1],'linewidth',2)
        plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'ro','markersize',12)
    end
    if ~isempty(Trace)
    nexttile([1 1])
    plot(Trace)
    hold all
    plot([i i],[min(Trace) max(Trace)],'color','r')
    axis tight
    end
    title([num2str(i/frmrate) ' ms'])
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
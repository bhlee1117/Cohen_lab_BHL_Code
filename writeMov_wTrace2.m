function writeMov_wTrace2(filename,movie,stack,rate,frmrate,range,dmd_mask_sequence_rois,Trace)

% this function generate a movie that showing the video and plot together. 
% Here the input video is X x Y x T and the input trace is X x T.
% If you want to show real time trace, check out writeMov_wTrace
% 2024.11.14 Byung Hun Lee

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
ylimval=[min(Trace(:)) max(Trace(:))];
%Trace=Trace(stack);
try
for i=1:length(stack)
    clf;
    tiledlayout(4,1)
    ax2=nexttile([3 1]);
    imshow2(movie(:,:,stack(i)),range)
    colormap('turbo')
    colorbar
    axis tight off equal
    hold on
    if ~isempty(dmd_mask_sequence_rois)
        %plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'color',[0 0.6 1],'linewidth',2)
        plot(dmd_mask_sequence_rois(:,1),dmd_mask_sequence_rois(:,2),'ro','markersize',12)
    end
    if ~isempty(Trace)
    ax1=nexttile([1 1]);
    plot(Trace(i,:))
    axis tight
    end
    ylim(ylimval)
    linkaxes([ax1 ax2],'x')
    %title([num2str(stack(i)/frmrate) ' ms'])
    %title([num2str(i*frmrate) '\mum'],'HorizontalAlignment','left')
    title([counting_string(i), ' lap'],'HorizontalAlignment','left')
    pause(0.005) %Pause and grab frame

    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    pause(0.05)
end
end
close(myVideo);
close(figure(1));
end
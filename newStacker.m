function S = newStacker(roi, cls, win, ev, H, W)
% NEWSTACKER  Initialise one (neuron, class) STA accumulator for S3/S4.
%   S = newStacker(roi, cls, win, ev, H, W)
%   roi : ROI/neuron index      cls : 1=SS, 2=CS
%   win : [pre post] in frames  ev  : event (trigger) frames for this modality
%   H,W : movie frame size
% Holds a running average sum plus a raw-window buffer that is flushed to
% numbered .bin parts by accumStack/flushStacker.
S.roi=roi; S.cls=cls; S.win=win; S.ev=ev(:)';
S.H=H; S.W=W;
S.sum=zeros(H,W,win(1)+win(2)+1);   % running sum of processed event windows (for the average)
S.cnt=0;                            % events added to the average
S.buf=[];                           % [Npx x winLen x nBuffered] raw windows not yet written
S.centers=[];                       % global center frames of buffered windows (in order)
S.centerFrames=[];                  % finalized ordered center frames (written parts)
S.binFiles={};                      % written .bin part filenames
S.evPerFile=[];                     % # windows in each written part
S.part=0;
end

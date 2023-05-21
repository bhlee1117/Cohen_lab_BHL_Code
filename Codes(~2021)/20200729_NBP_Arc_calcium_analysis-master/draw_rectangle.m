function draw_rectangle(crop_pos,linewidth,color,a)
if nargin<4
    a=gca;
end
line(a,[crop_pos(1,1) crop_pos(1,1)+crop_pos(1,3)],[crop_pos(1,2) crop_pos(1,2)],'linewidth',linewidth,'color',color)
line(a,[crop_pos(1,1) crop_pos(1,1)],[crop_pos(1,2) crop_pos(1,2)+crop_pos(1,4)],'linewidth',linewidth,'color',color)
line(a,[crop_pos(1,1) crop_pos(1,1)+crop_pos(1,3)],[crop_pos(1,2)+crop_pos(1,4) crop_pos(1,2)+crop_pos(1,4)],'linewidth',linewidth,'color',color)
line(a,[crop_pos(1,1)+crop_pos(1,3) crop_pos(1,1)+crop_pos(1,3)],[crop_pos(1,2) crop_pos(1,2)+crop_pos(1,4)],'linewidth',linewidth,'color',color)
end
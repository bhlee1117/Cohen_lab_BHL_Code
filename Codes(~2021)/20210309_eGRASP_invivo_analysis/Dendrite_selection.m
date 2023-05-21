function D_coor=Dendrite_selection(Im,c_image_xy,pix_xyz,method,D_coor)

i=[];
switch method
    case 'new'
        D_coor={[]};
    case 'add'

end
while isempty(i)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
     'Color',[1 1 1],'Renderer','painters','position',[1100 100 1200 1000]);
xz_rat=pix_xyz(1,2)/pix_xyz(1,1);

[sz1 sz2 sz3 sz4]=size(Im);
handles.axes1 = axes('Units','pixels','Position',[50 350 600 600]);  
imagesc(c_image_xy)
    hold all
    axis equal tight off
    cmap=distinguishable_colors(length(D_coor)-1);
    for coor=2:length(D_coor)
        plot(D_coor{coor}(:,1),D_coor{coor}(:,2),'marker','.','markersize',20,'color',cmap(coor-1,:),'linewidth',2)
    end
    
r=drawpolyline('color','r');
theta=rad2deg(atan((r.Position(end,2)-r.Position(1,2))/(r.Position(end,1)-r.Position(1,1))));
if theta<0
rot_im=imrotate3(squeeze(Im(:,:,:,3)),-abs(theta),[0 0 1],'loose');
else
rot_im=imrotate3(squeeze(Im(:,:,:,3)),abs(theta),[0 0 1],'loose');
end
[rotsz1 rotsz2 rotsz3]=size(rot_im);
clear r_rot inv_rot
for p=1:size(r.Position,1)
r_rot(p,:)=([cosd(-theta) -sind(-theta); sind(-theta) cosd(-theta)]*(r.Position(p,:)-[round(sz1/2) round(sz2/2)])')'+[round(rotsz1/2) round(rotsz2/2)];
end
hold all
r_rot
plot(r_rot(:,1),r_rot(:,2),'c')
%handles.axes2 = axes('Units','pixels','Position',[750 350 600*sz3*xz_rat/sz2 600]); 
%Seg_im=Im(:,min(round(r.Position(:,1))):max(round(r.Position(:,1))),:,3);
% rot_im=imrotate3(squeeze(Im(:,:,:,3)),-theta,[0 0 1]);
% Seg_im=rot_im(:,min(round(r_rot(:,1))):max(round(r_rot(:,1))),:);
% imagesc(fliplr(squeeze(max(Seg_im,[],2))))
% colormap('gray')
%  axis tight off
% h = images.roi.Polyline(gca,'Position',[repmat(10,size(r.Position,1),1) r.Position(:,2)],'color','r');

handles.axes3 = axes('Units','pixels','Position',[50 50 900 600*sz3*xz_rat/sz1]);  
%Seg_im=Im(min(round(r.Position(:,2))):max(round(r.Position(:,2))),:,:,3);

g=30;
try
Seg_im=rot_im(min(round(r_rot(:,2)))-g:max(round(r_rot(:,2)))+g,:,:);
catch 
    while (min(round(r_rot(:,2)))-g<1 | max(round(r_rot(:,2)))+g>sz1)
        g=g-1;
    end
    Seg_im=rot_im(min(round(r_rot(:,2)))-g:max(round(r_rot(:,2)))+g,:,:);
end

imagesc(imrotate(squeeze(max(Seg_im,[],1)),90))
colormap('gray')
 axis tight off
h2 = images.roi.Polyline(gca,'Position',[ r_rot(:,1) repmat(10,size(r_rot,1),1)],'color','r');

% handles.axes4 = axes('Units','pixels','Position',[750 310 600*sz3*xz_rat/sz2 30]);  
% imagesc(fliplr(reshape(jet(100),1,100,3)))
%  axis tight off
% 
handles.axes5 = axes('Units','pixels','Position',[970 50 30 600*sz3*xz_rat/sz1]);  
imagesc(flipud(reshape(jet(100),100,1,3)))
axis tight off

i=input(['Done = press enter to choose more or press other key to end\n']);
% if length(unique(h.Position(:,1)))>length(unique(h2.Position(:,2))) %use h
%     D_coor=[D_coor [r.Position sz3-interp1(h.Position(:,2),h.Position(:,1),r.Position(:,2),'pchip')]];
% else % use h2
for p=1:size(h2.Position,1)
%r_rot(p,:)=([cosd(-theta) -sind(-theta); sind(-theta) cosd(-theta)]*(r.Position(p,:)-[round(sz1/2) round(sz2/2)])')'+[round(rotsz1/2) round(rotsz2/2)];
inv_rot(p,:)=([cosd(theta) -sind(theta); sind(theta) cosd(theta)]*([h2.Position(p,1) mean(r_rot(:,2))]-round([rotsz1/2 rotsz2/2]))')'+[round(sz1/2) round(sz2/2)];
end



    D_coor=[D_coor, [r.Position sz3-interp1(h2.Position(:,1),h2.Position(:,2),r_rot(:,1))]];

% end

close all
end

switch method
    case 'new'
D_coor=D_coor(2:end);
    case 'add'
        if isempty(D_coor{1})
 D_coor=D_coor(2:end);  
        end
end


cmap=distinguishable_colors(length(D_coor));
for i=1:length(D_coor)
plot3(D_coor{i}(:,1)/xz_rat,D_coor{i}(:,2)/xz_rat,D_coor{i}(:,3),'marker','.','markersize',20,'color',cmap(i,:),'linewidth',2)
hold all
end
grid on
axis tight
  end




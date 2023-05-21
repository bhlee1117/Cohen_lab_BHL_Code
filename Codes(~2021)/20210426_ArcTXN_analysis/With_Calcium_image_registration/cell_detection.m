function cell_list=cell_detection(im_aligned,dist_thr,method)
N=2;

switch method
    case 'Find'
for i=1:size(im_aligned,2)
max_sub(:,:,i)=max(im_aligned{i},[],3);
[det{i}, radii]=Cell_segment_20200806(max_sub(:,:,i));
num_det(i,1)=size(det{i},1);
end
[max_cell_numb ref]=max(num_det(:,1));

zstack_circle=[];
h=waitbar(0,'Wait for XZ detection');
clear subt_crop

for i=1:size(im_aligned{ref},2)
    reslice(:,:,i)=im_aligned{ref}(:,i,:);
    [centers, radii]=Cell_segment_circle_zslice_20200806(reslice(:,:,i));
    if ~isempty(centers)
        centers=mod_cell_cent(centers,3);
        col=zeros(size(centers,1),1)+i;
zstack_circle(size(zstack_circle,1)+1:size(zstack_circle,1)+size(centers,1),1:3)=[centers(:,2) col centers(:,1)];
    end
    waitbar((i)/size(im_aligned{ref},2))
end
close(h)

for i=ref
    clear dist cell_list centers
    max_sub(:,:,i)=max(im_aligned{ref}(:,:,1:end),[],3);    
    g=1;
    h=waitbar(0,'Matching Cell position of Z and XY');
    zstack=size(im_aligned{ref},3);
    for n=1:N
max_sub_{1,n}(:,:,i)=max(im_aligned{ref}(:,:,(n-1)*floor(zstack/N)+1:n*floor(zstack/N)),[],3);
    [centers{1,n}, radii]=Cell_segment_20200806(max_sub_{1,n}(:,:,i));
    centers{1,n}=mod_cell_cent(centers{1,n},5);

    for j=1:size(zstack_circle,1)
        for k=1:size(centers{1,n},1)
            dist(k,j)=distance_BH([centers{1,n}(k,2) centers{1,n}(k,1)],zstack_circle(j,1:2));
        end
    end
        for k=1:size(centers{1,n},1)
        [min_dist min_dist_arg]=min(dist(k,:));
        if min_dist<dist_thr && zstack_circle(min_dist_arg,3)<n*round(zstack/N) && zstack_circle(min_dist_arg,3)>round(zstack/N)*(n-1)
            cell_list(g,1:3)=[centers{1,n}(k,1) centers{1,n}(k,2) zstack_circle(min_dist_arg,3)];
            g=g+1;
        end
        end
        
        waitbar((n)/N)
    end
    close(h)
end
    case 'Add'
        [fnm pth]=uigetfile('*.xls','Select Cell list file exported from Imaris');
         cell_list=readmatrix([pth fnm]);
         cell_list=cell_list(:,1:3);
end
    cell_list=mod_cell_cent(cell_list,5);
cmap=distinguishable_colors(size(im_aligned,2));
merge_image(im_aligned,cmap,[])
    hold all
    axis equal tight
    plot3(cell_list(:,1),cell_list(:,2),cell_list(:,3),'r.','markersize',20)
    x=1;
  title({'Press left click to add new cell, right click to remove existing cell'; 'Press enter to exit'},'fontweight','bold');
[x y button]=ginput(1);
close all
  while ~isempty(x)
    if button==1
      reslice(:,:,round(y))=im_aligned{1}(round(y),:,:);
      imagesc(imrotate(reslice(:,:,round(y)),90),[0 20])
      %imagesc(reshape(reslice(round(y),:,:),size(reslice,2),size(reslice,3)))
      colormap('gray')
      axis equal tight
      hold all
      draw_rectangle([x-20 1 40 size(im_aligned{1},3)],2,[1 0 0])
      [x z_manual]=ginput(1);
   %   z_manual=input('Input the z centroid position\n');
      close all
      cell_list(size(cell_list,1)+1,1:3)=[x y z_manual];
         elseif button==3
        for i=1:size(cell_list,1)
            d(i,1)=distance_BH(cell_list(i,1:2),[x y]);
        end
        [mind argd]=min(d);
        try
        cell_list(argd,:)=[];
        catch
        end
    end
    
merge_image(im_aligned,cmap,[])
      axis equal tight
      colormap('gray')
      hold all
      plot(cell_list(:,1),cell_list(:,2),'r.','markersize',20)
        title({'Press left click to add new cell, right click to remove existing cell'; 'Press enter to exit'},'fontweight','bold');
[x y button]=ginput(1);
      close all
 
  end
end
function [centers, radii]=Cell_segment_20200806(I)
if nargin<4
Cell_size=200; %pixels
kernel=2;
vac_ratio=0.6;
end

imfilt=double(imgaussfilt(I, kernel));
imfilt=imfilt*255/max(reshape(imfilt,size(imfilt,1)*size(imfilt,2),1));
imfilt=uint8(imfilt);

imfilt=2^8-1-imfilt;
%circleFinder(imfilt)W
detectCircles = @(x) imfindcircles(x,[14 17], ...
	'Sensitivity',0.9000, ...
	'EdgeThreshold',0.09, ...
	'Method','TwoStage', ...
	'ObjectPolarity','Dark');
[centers, radii, metric] = detectCircles(imfilt);

end
function [centers, radii]=Cell_segment_circle_zslice_20200806(I)
%% Control Tower
if nargin<4
Cell_size=200; %pixels
kernel=3;
vac_ratio=0.6;
end

imfilt=double(imgaussfilt(I, kernel));
imfilt=imfilt*255/max(reshape(imfilt,size(imfilt,1)*size(imfilt,2),1));
imfilt=uint8(imfilt);

imfilt=2^8-1-imfilt;
%circleFinder(imfilt)

detectCircles = @(x) imfindcircles(x,[17 20], ...
	'Sensitivity',0.9000, ...
	'EdgeThreshold',0.09, ...
	'Method','TwoStage', ...
	'ObjectPolarity','Dark');
[centers, radii, metric] = detectCircles(imfilt);

end

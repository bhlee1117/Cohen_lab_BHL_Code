%% Load the image and set the output path
clear
[fnm,pth]=uigetfile('*.tif','Select the timelapse images','Multiselect','on');
zstack =34; 
time=31;
dist_thr=3;
%% Make the cell index
    clear dist cell_list centers
    for i=1:zstack*5
        im(:,:,i)=imread([pth fnm],i);
    end
    %%
zstack_circle=[];
h=waitbar(0,'Detection in reslice imaging');
clear subt_crop
for i=1:zstack
    subt_crop(:,:,i)=imread([pth fnm],i);
end

for i=1:size(subt_crop,2)
    reslice(:,:,i)=subt_crop(:,i,:);
    
    [centers, radii]=Cell_segment_circle_zslice_512(reslice(:,:,i));
    if ~isempty(centers)
        centers=mod_cell_cent(centers,3);
        col=zeros(size(centers,1),1)+i;
zstack_circle(size(zstack_circle,1)+1:size(zstack_circle,1)+size(centers,1),1:3)=[centers(:,2) col centers(:,1)];

    end
    waitbar((i)/size(subt_crop,2))
end
close(h)
   max_sub(:,:,1)=max(subt_crop(:,:,1:zstack),[],3);    
   %%
    g=1; clear centers
    h=waitbar(0,'Matching Cell position of Z and XY');
    for n=1:2
max_sub_{1,n}(:,:,1)=max(im(:,:,(n-1)*floor(zstack/2)+1:n*floor(zstack/2)),[],3);
%max_sub_2(:,:,i)=max(subt_crop(:,:,41:81),[],3);
    [centers{1,n}, radii]=Cell_segment_circle_512(max_sub_{1,n}(:,:,1));
    centers{1,n}=mod_cell_cent(centers{1,n},5);
%     [centers2, radii]=Cell_segment_circle(max_sub_2(:,:,i));
%     centers2=mod_cell_cent(centers2,3);

    for j=1:size(zstack_circle,1)
        for k=1:size(centers{1,n},1)
            dist(k,j)=distance_BH([centers{1,n}(k,2) centers{1,n}(k,1)],zstack_circle(j,1:2));
        end
    end
        for k=1:size(centers{1,n},1)
        [min_dist min_dist_arg]=min(dist(k,:));
        if min_dist<dist_thr && zstack_circle(min_dist_arg,3)<n*round(zstack/2) && zstack_circle(min_dist_arg,3)>round(zstack/2)*(n-1)
            cell_list(g,1:3)=[centers{1,n}(k,1) centers{1,n}(k,2) zstack_circle(min_dist_arg,3)];
            g=g+1;
        end
        end
        
        waitbar((n)/2)
    end
    close(h)

    cell_list=mod_cell_cent(cell_list,7);
    %centers=Cell_segment2(max_sub(:,:,i),200,3,0.5);
    imagesc(max_sub(:,:,1))
    colormap('gray')
    hold all
    plot3(cell_list(:,1),cell_list(:,2),cell_list(:,3),'ro')
   % plot(centers1(:,1),centers1(:,2),'bo')
% mod_Cell_cent=mod_cell_cent(centers,15);
%     Cell_cent{1,i}=mod_Cell_cent;
%  

sw=menu('OK?','OK','More');
  while sw==2
      
      figure(2)
      imagesc(max_sub(:,:,ref))
      colormap('gray')
      hold all
      plot(cell_list(:,1),cell_list(:,2),'ro')
      [x y]=ginput(1);
      z=input('Input the z centroid position\n');
      cell_list(size(cell_list,1)+1,1:3)=[x y z];
      close(figure(2))
      
      figure(2)
      imagesc(max_sub(:,:,ref))
      colormap('gray')
      hold all
      plot(cell_list(:,1),cell_list(:,2),'ro')
      sw=menu('OK?','OK','More');
  end
Cells.list=cell_list;
Cells.filename=fnm;
Cells.max_image=max_sub;
Cells.zstack=zstack;
Cells.time=time;
save([pth 'Cell_list_' fnm '.mat'],'Cells');
%%
aviobj = VideoWriter([pth,fnm,'_seg.avi']);
open(aviobj);
rod_z_cell_list=cell_list;
rod_z_cell_list(:,3)=round(cell_list(:,3));
for i=1:zstack
    figure(1)
    imagesc(subt_crop(:,:,i))
    colormap('gray')
    axis equal
    hold all
    [row col]=find(rod_z_cell_list(:,3)==i);
    for j=1:size(row,1)
    plot(rod_z_cell_list(row(j,1),1),rod_z_cell_list(row(j,1),2),'r.','markersize',15)
    end
    
    F=figure(1);
writeVideo(aviobj,getframe(F));
end
close(figure(1));
 close(aviobj);
 %%
 aviobj = VideoWriter([pth,fnm,'_seg_resl.avi']);
open(aviobj);
rod_z_cell_list=cell_list;
rod_z_cell_list(:,1)=round(cell_list(:,1));

for i=1:size(subt_crop,2)
    figure(1)
    imagesc(reshape(subt_crop(:,i,:),size(subt_crop,1),size(subt_crop,3)))
    colormap('gray')
    axis equal
    hold all
    [row col]=find(rod_z_cell_list(:,1)==i);
    for j=1:size(row,1)
    plot(rod_z_cell_list(row(j,1),3),rod_z_cell_list(row(j,1),2),'r.','markersize',15)
    end
    
    
    F=figure(1);
writeVideo(aviobj,getframe(F));
if ~mod(i,100)
        close(figure(1))
    end
end
close(figure(1));
 close(aviobj);
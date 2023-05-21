%%
function cell_list=G_cell_find(im,ratio)
zstack_circle=[];
res_im=imresize(im,ratio);
for i=1:size(res_im,2)
    
    reslice(:,:,i)=res_im(:,i,:);
    [centers, radii]=Cell_segment_circle_zslice_mem_act(reslice(:,:,i));
    if ~isempty(centers)
        centers=mod_cell_cent(centers,3);
        col=zeros(size(centers,1),1)+i;
zstack_circle(size(zstack_circle,1)+1:size(zstack_circle,1)+size(centers,1),1:3)=[centers(:,2) col centers(:,1)];

    end
end

    clear dist cell_list centers
    dist_thr=7;
    g=1;
   n=1; N=1;
max_res_im=max(res_im,[],3);
    [centers, radii]=Cell_segment_circle_mem_act(max_res_im);
    centers=mod_cell_cent(centers,5);

    for j=1:size(zstack_circle,1)
        for k=1:size(centers,1)
            dist(k,j)=distance_BH([centers(k,2) centers(k,1)],zstack_circle(j,1:2));
        end
    end
        for k=1:size(centers,1)
        [min_dist min_dist_arg]=min(dist(k,:));
        if min_dist<dist_thr && zstack_circle(min_dist_arg,3)<n*round(size(res_im,3)/N) && zstack_circle(min_dist_arg,3)>round(size(res_im,3)/N)*(n-1)
            cell_list(g,1:3)=[centers(k,1) centers(k,2) zstack_circle(min_dist_arg,3)];
            g=g+1;
        end
        end
   

    cell_list=mod_cell_cent(cell_list,7);
    %centers=Cell_segment2(max_sub(:,:,i),200,3,0.5);
     imagesc(max_res_im)
    colormap('gray')
    hold all
    plot3(cell_list(:,1),cell_list(:,2),cell_list(:,3),'ro')

sw=menu('OK?','OK','More');
  while sw==2
      
      figure(2)
      imagesc(max_res_im)
      colormap('gray')
      hold all
      plot(cell_list(:,1),cell_list(:,2),'ro')
      [x y]=ginput(1);
      z=input('Input the z centroid position\n');
      cell_list(size(cell_list,1)+1,1:3)=[x y z];
      close(figure(2))
      
      figure(2)
      imagesc(max_res_im)
      colormap('gray')
      hold all
      plot(cell_list(:,1),cell_list(:,2),'ro')
      sw=menu('OK?','OK','More');
  end

close(figure(2))
close(figure(1))
cell_list(:,1:2)=cell_list(:,1:2)/ratio;
end

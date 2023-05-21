function cell_list=cell_detection_manual(im,cell_list,rng)
N=2;
f=figure('units','normalized','outerposition',[0 0 1 1]);
imshow2(im,rng);
    hold all
    axis equal tight
    colormap('gray')
    plot(cell_list(:,1),cell_list(:,2),'r.','markersize',15)
    x=1;
  title({'Press left click to add new cell, right click to remove existing cell'; 'Press enter to exit'},'fontweight','bold');
[x y button]=ginput(1);
close(f)
  while ~isempty(x)
    if button==1
   %   z_manual=input('Input the z centroid position\n');
      %close all
      cell_list(size(cell_list,1)+1,1:2)=[x y];
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
    f=figure('units','normalized','outerposition',[0 0 1 1]);
    imshow2(im,rng);
    hold all
    axis equal tight
    colormap('gray')
    plot(cell_list(:,1),cell_list(:,2),'r.','markersize',15)
        title({'Press left click to add new cell, right click to remove existing cell'; 'Press enter to exit'},'fontweight','bold');
[x y button]=ginput(1);
      close(f)
 
  end
end

function distt=distance_BH(X,Y)
for i=1:size(X,1)
    distt(i,1)=0;
    for j=1:size(X,2)
    distt(i,1)=distt(i,1)+(X(i,j)-Y(i,j))^2;
    end
    distt(i,1)=sqrt(distt(i,1));
end
end

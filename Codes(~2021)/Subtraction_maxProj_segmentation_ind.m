%%%  TXN analysis protocol : Subtraction & Max projection & Cell Segmentation 
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code shows cell image and let user classify the TS cell and no TS
% cell.
% The analysis is proceed as following
% 1. Load the registrated images
% 2. Do subtraction in each z stack
% 3. Max projection
% 4. Cell Segment with Circle finder.

% INPUTS 

% The folder with images (3 stack, Reg, Green, Subtracted).
% 

% OUTPUTS


% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.

%% Parameter setting
clear
[fnm,pth]=uigetfile('*.tif','Select the registerated HC, Train and retrieval images','Multiselect','on');
channel=2;
timeline={'HC','Training','Retrieval'};
dist_thr=7;
N=2;
% target_RGB_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop\t01_z01_cell_num0019.tif';
% target_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop_sub\t01_z01_cell_num0019.tif';

%%
mkdir([pth,'crop'])
for i=1:length(fnm)
    
    
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);    
zstack=numstack/channel;

for j=1:zstack
    im(:,:,1)=imread([pth char(fnm{1,i})],j);
    imG=imread([pth char(fnm{1,i})],j+zstack);
    
    pseudo=imG>2500;
    imG(pseudo)=median(reshape(imG,size(imG,1)*size(imG,2),1));
    im(:,:,2)=imG;
    
subt_crop(:,:,j)=subtraction_image_ftn(im);
if j==1
    imwrite(subt_crop(:,:,j),[pth '\crop\' 'sub_' char(fnm{1,i})],'Compression','none');
else
imwrite(subt_crop(:,:,j),[pth '\crop\' 'sub_' char(fnm{1,i})],'Writemode','Append','Compression','none');
end
end

max_sub(:,:,i)=max(subt_crop,[],3);
[det{i}, radii]=Cell_segment_circle(max_sub(:,:,i));
num_det(i,1)=size(det{i},1);

end
%%

h=waitbar(0,'Processing');
clear subt_crop
for j=3%:length(fnm)
    zstack_circle{j}=[];
for i=1:zstack
    subt_crop{j}(:,:,i)=imread([pth '\crop\' 'sub_' char(fnm{1,j})],i);
end

for i=1:size(subt_crop{j},2)
    reslice{j}(:,:,i)=subt_crop{j}(:,i,:);
    
    [centers, radii]=Cell_segment_circle_zslice(reslice{j}(:,:,i));
    if ~isempty(centers)
        centers=mod_cell_cent(centers,3);
        col=zeros(size(centers,1),1)+i;
zstack_circle{j}(size(zstack_circle{j},1)+1:size(zstack_circle{j},1)+size(centers,1),1:3)=[centers(:,2) col centers(:,1)];

    end
    
end
waitbar((j)/length(fnm))
end
close(h)
%% Cell Segmentation in each image and match
for i=3%:length(fnm)
    clear dist cell_list centers
    max_sub(:,:,i)=max(subt_crop{i}(:,:,1:zstack),[],3);    
    g=1;
    h=waitbar(0,'Processing');
    clear centers
    for n=1:N
max_sub_{1,n}(:,:,i)=max(subt_crop{i}(:,:,(n-1)*floor(zstack/N)+1:n*floor(zstack/N)),[],3);
%max_sub_2(:,:,i)=max(subt_crop(:,:,41:81),[],3);
    [centers{1,n}, radii]=Cell_segment_circle(max_sub_{1,n}(:,:,i));
    centers{1,n}=mod_cell_cent(centers{1,n},5);
%     [centers2, radii]=Cell_segment_circle(max_sub_2(:,:,i));
%     centers2=mod_cell_cent(centers2,3);

    for j=1:size(zstack_circle{i},1)
        for k=1:size(centers{1,n},1)
            dist(k,j)=distance_BH([centers{1,n}(k,2) centers{1,n}(k,1)],zstack_circle{i}(j,1:2));
        end
    end
        for k=1:size(centers{1,n},1)
        [min_dist min_dist_arg]=min(dist(k,:));
        if min_dist<dist_thr && zstack_circle{i}(min_dist_arg,3)<n*round(zstack/N) && zstack_circle{i}(min_dist_arg,3)>round(zstack/N)*(n-1)
            cell_list{i}(g,1:3)=[centers{1,n}(k,1) centers{1,n}(k,2) zstack_circle{i}(min_dist_arg,3)];
            g=g+1;
        end
        end
        
        waitbar((n)/N)
    end
    close(h)
    
%     
%     for j=1:size(zstack_circle,1)
%         for k=1:size(centers2,1)
%             dist(k,j)=distance_BH([centers2(k,2) centers2(k,1)],zstack_circle(j,1:2));
%         end
%     end
%         for k=1:size(centers2,1)
%         [min_dist min_dist_arg]=min(dist(k,:));
%         if min_dist<dist_thr && zstack_circle(min_dist_arg,3)>40
%             cell_list(g,1:3)=[centers2(k,1) centers2(k,2) zstack_circle(min_dist_arg,3)];
%             g=g+1;
%         end
%         end
%     
    cell_list{i}=mod_cell_cent(cell_list{i},7);
    %centers=Cell_segment2(max_sub(:,:,i),200,3,0.5);
    imagesc(max_sub(:,:,i))
    colormap('gray')
    hold all
    plot3(cell_list{i}(:,1),cell_list{i}(:,2),cell_list{i}(:,3),'ro')
   % plot(centers1(:,1),centers1(:,2),'bo')
% mod_Cell_cent=mod_cell_cent(centers,15);
%     Cell_cent{1,i}=mod_Cell_cent;
%  

sw=menu('OK?','OK','More');
  while sw==2
      
      figure(2)
      imagesc(max_sub(:,:,i))
      colormap('gray')
      hold all
      plot(cell_list{i}(:,1),cell_list{i}(:,2),'ro')
      [x y]=ginput(1);
      z=input('Input the z centroid position\n');
      cell_list{i}(size(cell_list{i},1)+1,1:3)=[x y z];
      close(figure(2))
      
      figure(2)
      imagesc(max_sub(:,:,i))
      colormap('gray')
      hold all
      plot(cell_list{i}(:,1),cell_list{i}(:,2),'ro')
      sw=menu('OK?','OK','More');
  end
  Cells.list=cell_list{i};
Cells.filename=char(fnm{1,i});
Cells.max_image=max_sub(:,:,i);
Cells.zstack=zstack;
save([pth 'Cell_list_' char(fnm{1,i}) '.mat'],'Cells');
end

%%
aviobj = VideoWriter([pth,char(fnm{1,2}),'_seg.avi']);
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
 aviobj = VideoWriter([pth,char(fnm{1,2}),'_seg_resl.avi']);
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
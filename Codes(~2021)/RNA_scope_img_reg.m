%% Functional Role : 1. Image realignment for shading correction and camera shift 2. Stitching the Grid image
% Wroten by Byung Hun Lee 16/01/2018
%% Parameter Setting

clear; clc; close all;
Pathname='H:\Image_data\WF\20191119_RNA_scope\20191119_CFC2-1_4\'; %Image path
GridN=4;
tot_pix=512;
overlap_pix=50;
zmin=1;
zmax=11;
ch=4;
save_dir=Pathname;

%% Max projection
g=1;
clear im
for i=0:GridN-1
    for j=0:GridN-1
            
Pathname2=['20191119_CFC2-1_4_MMStack_2-Pos_',num2str(j,'%03d'),'_',num2str(i,'%03d'),'.ome'];
for c=1:ch
for k=zmin:zmax
   im{c}(:,:,k)=imread([Pathname,Pathname2,'.tif'],k+(c-1)*(zmax-zmin+1)); 
end
maxim{c}=max(im{c},[],3);

imwrite(maxim{c},[Pathname,'MAX_',Pathname2,'.tif'],'Writemode','append')
end
baSicim(:,:,g)=maxim{1};
g=g+1;
%imwrite(max(imR,[],3),[Pathname,'MAX_',Pathname2,'.tif'],'WriteMode','append')

    end
end

[flatfield_G,darkfield_G] = BaSiC(baSicim,'darkfield','true');  
field_TIRF(:,:,1)=flatfield_G;
field_TIRF(:,:,2)=darkfield_G;
save([save_dir,'field.mat'],'field_TIRF')

%%
% load ([matrixpath, tformFile]);
% load([save_dir,'field.mat'])
clear map
sqr=tot_pix-overlap_pix;
h=waitbar(0,'Constructing the corrected map');
for c=1:ch
im_corr_G{c} = zeros(size(im{c}));
for i=0:GridN-1
    for j=0:GridN-1
        Pathname2=['20191119_CFC2-1_4_MMStack_2-Pos_',num2str(j,'%03d'),'_',num2str(i,'%03d'),'.ome'];
        name_f=['MAX_' Pathname2];
        
    
    %im_corr_G{c}(:,:,i*GridN+j+1) = (double(imread([Pathname,name_f,'.tif'],c))-darkfield_G)./flatfield_G;
        im_corr_G{c}(:,:,i*GridN+j+1) = (double(imread([Pathname,name_f,'.tif'],c)));
    
    map{c}(i*sqr+1:(i+1)*sqr,j*sqr+1:(j+1)*sqr)=fliplr(im_corr_G{c}(overlap_pix/2+1:overlap_pix/2+sqr,overlap_pix/2+1:overlap_pix/2+sqr,i*GridN+j+1));
    
    end
    end
end

close(h)

% map_G(j*412+1:j*412+412,i*412+1:i*412+412)=fliplr(maxG(31:442,31:442));
% map_R(j*412+1:j*412+412,i*412+1:i*412+412)=fliplr(maxR(31:442,31:442));
for c=1:ch    
imwrite(uint16(map{c}),[save_dir,'map.tif'],'WriteMode','append');
end


%% Load files
clear
[fnm pth]=uigetfile('*.tif','Multiselect','on');
ch=4;
%% Cell detection
if ischar(fnm)
    fnm={fnm};
end

for i=1:length(fnm)
    clear im im_tmp maxIm RNAscope
    iminfo=imfinfo([pth fnm{1,i}]);
    ctrd{i}=[];
    for c=1:ch
    for j=1:size(iminfo,1)
    im_tmp=imread([pth fnm{1,i}],j);  
    
    im{c}(:,:,j)=im_tmp(:,:,c);
    end
    maxIm{i,c}=max(im{c},[],3);
    RNAscope.max(:,:,c)=maxIm{i,c};
    end
    
    [centers, radii]=Cell_segment_circle_RNA_scope_20191120(maxIm{i,1});
    ctrd{i}=[ctrd{i}; centers];
    ctrd_md{i}=mod_cell_cent(ctrd{i},40);
    [tm tm_rank]=sort(ctrd_md{i}(:,2),'descend');
    RNAscope.list=ctrd_md{i}(tm_rank,:);
        
save([pth fnm{1,i} '.mat'],'RNAscope')
end

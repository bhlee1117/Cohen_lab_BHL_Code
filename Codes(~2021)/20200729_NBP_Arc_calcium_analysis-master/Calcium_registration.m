function Calcium_registration(pth,ch,ref,step,save_tif)
  % Select the .tif file folder taken from the olympus Two photon microscope. 
%%
for p=1:length(pth)
if ch~=1 && ref==2
    ap=1;
else
    ap=2;
end
clear order_img order
    img=imageDatastore(pth{1,p});
    sp=split(img.Files,'T');
    sp=split(sp(:,2),'.');
    for z=1:size(sp,1)
    order_img(z,1)=str2num(sp{z,1});
    end
    [as order]=sort(order_img,'ascend');
    %step=floor(size(im{1},3)/slice);
    if isempty(step)
    step=size(sp,1)/ch;
    end
    slice=floor(size(order_img,1)/step/ch);
    for j=1:slice
    switch ch
        case 1
    im{1}=[]; im{2}=[];
   
    for z=1:step
        im{1}(:,:,z)= (char(img.Files(order((step*(j-1)+z),1),:)));
    end
    Y = single(im{1});                 % convert to single precision 
    T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
    tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
    imwrite(uint16(mean(M1,3)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
        
        case 2
    im{1}=[]; im{2}=[];
     f = waitbar(0,'Loading images  ');
    for z=1:step
        waitbar(z/step,f,['Loading images  ' num2str(z) '/' num2str(step)]);
        im{1}(:,:,z)=imread(char(img.Files(order(2*(step*(j-1)+z)-1,1),:)));
        im{2}(:,:,z)=imread(char(img.Files(order(2*(step*(j-1)+z),1),:)));
    end
    close(f)
    Y = single(im{ref});                 % convert to single precision 
    T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
    tic; [Y,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
    MR(:,:,j)=mean(Y,3);
    MG(:,:,j) =mean(apply_shifts(im{ap},shifts1,options_rigid),3);
    end
    end
     
     save([pth{1,p}(1,1:end-10) '_calcium_image.mat'],'MR','MG','Y','-v7.3');
     if save_tif 
         for i=1:size(Y,3)
             
             imwrite(uint8(Y(:,:,i)),[pth{1,p}(1,1:end-10) '_calcium_imag.tif'],'Writemode','Append')
         end
         imwrite(uint16(M1(:,:,j)),[pth{1,p}(1,1:end-10) 'Rmean.tif'],'Writemode','Append')
     imwrite(uint16(M2(:,:,j)),[pth{1,p}(1,1:end-10) 'Gmean.tif'],'Writemode','Append')
     end
end

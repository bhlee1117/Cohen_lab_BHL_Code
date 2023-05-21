[fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
mov = vm(fpath);
% cd(fpath)
%% calculate centroid

[Y,X] = ndgrid(1:mov.rows,1:mov.cols);
mov_data = double(mov.data);
x_c = sum(mov_data .* X,[1 2])./sum(mov_data,[1 2]);
y_c = sum(mov_data .* Y,[1 2])./sum(mov_data,[1 2]);

figure;
imshow(mov.mean,[],'initialmagnification','fit')
hold on
scatter(x_c,y_c)
%%
figure;
for ii=1:100:length(x_c)
    clf
    imshow(mov.mean,[],'initialmagnification','fit')
    hold on
    scatter(x_c(ii),y_c(ii))
    title(sprintf('%g',ii))
    pause(.1)
end
%%
mov = double(mov.data);
mov([1 2 end],:,:)=0;
mov = vm(mov);
mov2 = mov.blur(3);
mov2 = mov2.data;
% figure;moviesc(mov)

[im_y,im_x] =ndgrid(1:mov.size(1),1:mov.size(2));
[spt_y_interp,spt_x_interp] =meshgrid(linspace(1,mov.size(1),500),linspace(1,mov.size(2),500));
spt_x = zeros(size(mov2,3),1);spt_y = zeros(size(mov2,3),1);
% figure
mov2 = gpuArray(mov2);
im_x = gpuArray(im_x); im_y = gpuArray(im_y);
spt_x_interp = gpuArray(spt_x_interp); spt_y_interp = gpuArray(spt_y_interp);
spt_x = gpuArray(spt_x); spt_y = gpuArray(spt_y);
for ii=1:mov.size(3)
%     clf
    im_interp = interp2(im_x,im_y,squeeze(mov2(:,:,ii)),spt_x_interp,spt_y_interp,'cubic');
    [~,spt_ind] = max(im_interp(:));
%     surf(im_interp,'edgecolor','none')
%     pause(.1)
%     [spt_y_ind,spt_x_ind] = ind2sub(size(spt_x_interp),spt_ind);
    spt_y(ii) = spt_y_interp(spt_ind);
    spt_x(ii) = spt_x_interp(spt_ind);
end

spt_x = gather(spt_x);
spt_y = gather(spt_y);
% spt_x = squeeze(mean(spt_x.*mov.data./sum(mov.data,[1 2])*mov.size(1)*mov.size(2),[1 2]));
% spt_y = squeeze(mean(spt_y.*mov.data./sum(mov.data,[1 2])*mov.size(1)*mov.size(2),[1 2]));


% [~,spt_ind] = max(mov(:,:));
% [spt_y,spt_x] = ind2sub(mov.size(1:2),spt_ind);



figure;
imshow(mov.mean,[],'initialmagnification','fit');
hold on
scatter(spt_x,spt_y,'o')
%%
figure(10);clf%;moviesc(mov)
subplot(3,2,1)
plot(spt_x)
title('spot center x')
subplot(3,2,2)
plot(spt_y)
title('spot center y')
subplot(3,2,3)
plot(squeeze(mean(mov.data,[1 2])))
subplot(3,2,4)
plot(squeeze(mean(mov.data,[1 2])))
subplot(3,2,5)
% histogram2(spt_x,spt_y,50,'DisplayStyle','tile','ShowEmptyBins','on');
% colorbar
nGrp = 100;
gscatter(spt_x,spt_y,round(linspace(1,nGrp,length(spt_x))),colormap(jet(nGrp)));legend off;daspect([.2 1 1])
c = colorbar;
c.Label.String = 'Frame No.';
c.TickLabels = linspace(0,length(spt_x),5);
c.Ticks = linspace(0,1,5);
% subplot(3,2,6)

% axis equal
% imshow(squeeze(mov.data(:,:,i)),[],'initialmagnification','fit')
%     hold on
%     plot(spt_x(i),spt_y(i),'o')
%%

for i = 1:10:length(spt_x)
    clf
    imshow(squeeze(mov.data(:,:,i)),[],'initialmagnification','fit')
    hold on
    plot(spt_x(i),spt_y(i),'o')
    
    title(sprintf('%g',i))
%     ylim([13 14]);xlim([27 29])
    pause(.1)
end




clear
[ fnm pth]=uigetfile('*.tif','Select the TXN data','Multiselect','on');
[fnm_ca pth_ca ]=uigetfile('*.tif','Select the Calcium image');
load('NF_map_.mat')
fnm_ca={fnm_ca};
fnm={fnm};
cls_pth='H:\Image_data\OlympusTPM\20190626_jRCaMP1a2\Classification';
sz=200;
video=0;
timescale=0.31; %Second
%%
clear dff
clear imexp
cell_plane=NF_map.cell;
sz=200;
Cell_size=[61 61 61];
ncell=size(find(NF_map.list(:,4)~=3 & NF_map.list(:,5)~=3)',2)+1;
nt=1;
t=1;

imginf=imfinfo([pth_ca fnm_ca{1,1}]); %%Generating ca image
numstack=size(imginf,1);
for i=1:numstack
    im_ca_crop(:,:,i)=imread([pth_ca char(fnm_ca)],i);
end

imexp=zeros(size(im_ca_crop,1)+2*floor(sz)/2,size(im_ca_crop,2)+2*floor(sz/2),numstack);
      for i=1:numstack
    imexp(sz/2+1:sz/2+size(im_ca_crop,1),sz/2+1:sz/2+size(im_ca_crop,2),i)=im_ca_crop(:,:,i);
      end
     cell_plane_shift=cell_plane+[sz/2 sz/2 Cell_size(1,3)/2];
    imginf=imfinfo([pth fnm{1,1}]); %%Generating ca image
numstac=size(imginf,1);
for i=1:numstac
    im_TXN_reg_crop(:,:,i)=imread([pth char(fnm)],i);
end

    im_exp=zeros(size(im_TXN_reg_crop,1)+2*floor(sz/2),size(im_TXN_reg_crop,2)+2*floor(sz/2),size(im_TXN_reg_crop,3)+2*floor(Cell_size(1,3)/2));
    im_exp(floor(sz/2)+1:floor(sz/2)+size(im_TXN_reg_crop,1),floor(sz/2)+1:floor(sz/2)+size(im_TXN_reg_crop,2),...
           floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im_TXN_reg_crop,3))=im_TXN_reg_crop;
       nt=1; t=1;
       
       for j=find(NF_map.list(:,4)~=3 & NF_map.list(:,5)~=3)'
            if NF_map.list(j,5)==1
       dff{1}(nt,:)=(NF_map.Calcium(j,:)-median(NF_map.Calcium(j,:)))/median(NF_map.Calcium(j,:));
       nt=nt+1;
            else
                 dff{2}(t,:)=(NF_map.Calcium(j,:)-median(NF_map.Calcium(j,:)))/median(NF_map.Calcium(j,:));
                 t=t+1;
            end
            
       end
       %%
    if video==1     
aviobj = VideoWriter([pth,'ca_trace.avi']);
open(aviobj);
st=1;
    else
        st=numstack;
    end
       % Figure generation
handles.fig = figure(1);
set(gcf, 'Position',  [50, 30, 1500, 950])
for i=st:numstack
    nt=1; t=1;
   
for j=find(NF_map.list(:,4)~=3 & NF_map.list(:,5)~=3)'
    
crop_im=im_exp(round(cell_plane_shift(j,2))-floor(Cell_size(1,2)/2):round(cell_plane_shift(j,2))+floor(Cell_size(1,2)/2),...
               round(cell_plane_shift(j,1))-floor(Cell_size(1,1)/2):round(cell_plane_shift(j,1))+floor(Cell_size(1,1)/2),...
               round(cell_plane_shift(j,3))-floor(Cell_size(1,3)/2):round(cell_plane_shift(j,3))+floor(Cell_size(1,3)/2));
           
    if NF_map.list(j,5)==1
       
handles.axes1 = axes('Units','pixels','Position',[50 870/ncell*nt+30 0.9*900/ncell 0.9*900/ncell]);    
    hold all
crop_im_ca(:,:,i)=imexp(round(cell_plane_shift(j,2))-sz/2:round(cell_plane_shift(j,2))+sz/2,round(cell_plane_shift(j,1))-sz/2:round(cell_plane_shift(j,1))+sz/2,i);
imshowpair(crop_im_ca(:,:,i),NF_map.mask(:,:,j),'ColorChannels','red-cyan')
axis equal tight off
handles.axes2 = axes('Units','pixels','Position',[100 870/ncell*nt+30 0.9*900/ncell 0.9*900/ncell]);
imagesc(max(crop_im,[],3))
colormap('gray')
hold all
axis equal tight off
axis off
nt=nt+1;
    else
        
handles.axes1 = axes('Units','pixels','Position',[50 870/ncell*(size(dff{1},1)+t)+45 0.9*870/ncell 0.9*870/ncell]);   
hold all
crop_im_ca(:,:,i)=imexp(round(cell_plane_shift(j,2))-sz/2:round(cell_plane_shift(j,2))+sz/2,round(cell_plane_shift(j,1))-sz/2:round(cell_plane_shift(j,1))+sz/2,i);
imshowpair(crop_im_ca(:,:,i),NF_map.mask(:,:,j),'ColorChannels','red-cyan')
axis equal tight off
handles.axes2 = axes('Units','pixels','Position',[100 870/ncell*(size(dff{1},1)+t)+45 0.9*870/ncell 0.9*870/ncell]);
imagesc(max(crop_im,[],3))
colormap('gray')
axis equal tight off
hold all
%axis off
t=t+1;     
    end
             

end
handles.axes3 = axes('Units','pixels','Position',[150 40 1300 900/ncell*size(dff{1},1)]);
handles.axes4 = axes('Units','pixels','Position',[150 900/ncell*size(dff{1},1)+40 1300 900/ncell*(size(dff{2},1)+1)]);
    
plot(handles.axes3,[timescale:timescale:timescale*i],dff{1}(:,1:i)+0.4*[1:1:size(dff{1},1)]')
set(handles.axes3,'YTick',[0.4:0.4:0.4*size(dff{1},1)],'YTickLabel',...
    [1:1:size(dff{1},1)]);
xlabel(handles.axes3,'Time (sec)')
hold all
xlim(handles.axes3,[timescale timescale*numstack])
ylim(handles.axes3,0.4*[0 size(dff{1},1)+1])
 
plot(handles.axes4,[timescale:timescale:timescale*i],dff{2}(:,1:i)+0.4*[1:1:size(dff{2},1)]')
hold all
set(handles.axes4,'YTick',[0.4:0.4:0.4*size(dff{2},1)],'YTickLabel',...
    [1:1:size(dff{2},1)]);
xlim(handles.axes4,[timescale timescale*numstack])
ylim(handles.axes4,0.4*[0.5 size(dff{2},1)+1.5])
F=figure(1);
if video==1
writeVideo(aviobj,getframe(F));
end
if mod(i,50)==0 && video==1
    close(figure(1))
    handles.fig = figure(1);
set(gcf, 'Position',  [50, 30, 1500, 950])

end

end
 close(figure(1));
 if video==1
 close(aviobj);
 end

%% OASIS
g=1;
for n=1:2
for i=1:size(dff{n},1)
[c_oasis, s_oasis,options] = deconvolveCa(dff{n}(i,:), 'ar1', 'constrained');
denoised{n}(i,:)=c_oasis;
deconv{n}(i,:)=s_oasis;
opt{n,i}=options;


plot([timescale:timescale:timescale*numstack],dff{n}(i,:)+g*0.5,'color','k')
hold all
plot([timescale:timescale:timescale*numstack],denoised{n}(i,:)+g*0.5,'color','r')
for j=1:size(dff{n},2)
    if deconv{n}(i,j)~=0
line([timescale*j timescale*j],[g*0.5 deconv{n}(i,j)+g*0.5],'color',[0 0 1])
hold all
    end
end
xlim([0 timescale*numstack])

g=g+1;
end
end
set(gca,'YTick',[0.5:0.5:0.5*(g-1)],'YTickLabel',...
    [1:1:g-1]);
NF_map.deconv=deconv;
NF_map.denois=denoised;
save([pth,'NF_map_.mat'],'NF_map');

%%
handles.fig = figure(1);
set(gcf, 'Position',  [50, 30, 1500, 950])
handle.axes1 = axes('Units','pixels','Position',[50 50 1500 450]);
imagesc(deconv{1},[0 0.2])
colormap('gray')
colorbar
handle.axes2= axes('Units','pixels','Position',[50 550 1500 450/size(dff{1},1)*size(dff{2},1)]);
imagesc(deconv{2},[0 0.2])
colormap('gray')
colorbar

handles.fig = figure(2);
set(gcf, 'Position',  [50, 30, 1500, 950])
handle.axes1 = axes('Units','pixels','Position',[50 50 1500 450]);
imagesc(dff{1},[0 0.2])
colormap('jet')
colorbar
handle.axes2= axes('Units','pixels','Position',[50 550 1500 450/size(dff{1},1)*size(dff{2},1)]);
imagesc(dff{2},[0 0.2])
colormap('jet')
colorbar

imwrite(uint16(deconv{1}*100),[pth 'non_txn.tif'],'Writemode','Append')
imwrite(uint16(denoised{1}*100),[pth 'non_txn.tif'],'Writemode','Append')
imwrite(uint16(dff{1}*100),[pth 'non_txn.tif'],'Writemode','Append')
imwrite(uint16(deconv{2}*100),[pth 'txn_deconv.tif'],'Writemode','Append')
imwrite(uint16(denoised{2}*100),[pth 'txn_deconv.tif'],'Writemode','Append')
imwrite(uint16(dff{2}*100),[pth 'txn_deconv.tif'],'Writemode','Append')
%%
T = 0.31;             % Sampling period       
Fs=1/T;
L = 1000;             % Length of signal
t = (0:L-1)*T; 

for i=1:size(dff,2)
    for j=1:size(dff{i},1)
        subplot(2,1,i)
FT_dff{i}(j,:)=fft(deconv{i}(j,:));
P2 = abs(FT_dff{i}(j,:)/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
% if i==1
% plot(f,P1,'color',[0.6 0.6 0.6]) 
% else
% plot(f,P1,'color','r') 
% alpha(0.1)
% end
plot(f,P1)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
hold all
    end
end

%%
for i=1:size(dff,2)
    for j=1:size(dff{i},1)
 
[aco_dff{i}(j,:) bin]=autocorr(deconv{i}(j,:),size(deconv{i},2)/2);
bin=bin*T;
if i==1
plot(bin,aco_dff{i}(j,:),'color',[0.6 0.6 0.6]) 
else
plot(bin,aco_dff{i}(j,:),'color','r') 
end
xlabel('Lag time (Sec)')
ylabel('Auto-correlation')
hold all
    end
end
%%
T = 0.31;             % Sampling period     
clear df sumsplike
df=[deconv{1};deconv{2}];

    for j=1:size(df,1)
sumspike(j,1)=sum(df(j,:));
    end

%%
clear df crscor arg
df=[deconv{1};deconv{2}];
gg=1;
g=1;
for i=find(sumspike>1)'
    g=1;
    for j=find(sumspike>1)'
[cor cor_bin]=crosscorr(df(i,:),df(j,:),'NumLags',499);
%[crscor(g,gg) arg(i,j)]=max(abs(cor));
[crscor(g,gg) ]=cor(1,500);
g=g+1;
    end
    gg=gg+1;
end
imagesc(crscor)
colormap('jet')
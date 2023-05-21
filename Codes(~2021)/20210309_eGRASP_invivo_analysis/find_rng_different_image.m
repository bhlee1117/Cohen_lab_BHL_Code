 function [candidatePos Y reg_angle shift ]=find_rng_different_image(D_coor,Im)
% Im(1,x) : Day -1 image, Im(2,x) : Day +2 image

m_image{1}=max(Im{1,1}(:,:,:,3),[],3); m_image{2}=max(Im{2,1}(:,:,:,3),[],3);
D_name={'Day -1','Day +2'};
[sz1 sz2]=size(m_image{1});
cmap= distinguishable_colors(size(D_coor,2));

figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
     'Color',[1 1 1],'Renderer','painters','position',[100 100 1000 500]);
 angle=[-5:0.2:10];
for a=1:length(angle)
    rot_im=imrotate(m_image{2},angle(a));
    [ match_data, ~ ] = matchPattern(rot_im, max(m_image{1}(round(sz1/2)-100:round(sz1/2)+100,round(sz2/2)-100:round(sz2/2)+100,:),[],3),...
                                     0.2,2); 
    [x xx]=max(match_data(:,3));
    corr_coef(a,:)=match_data(xx,:);
end
[m argm]= max(corr_coef(:,3)); reg_angle=angle(argm);
ch=length(size(Im{1,1}));
for c=1:ch
Y{1,c}=Im{1,1}(:,:,:,c); end
for z=1:size(Im{2,1},3)
    for c=1:ch
    Y{2,c}(:,:,z)=imrotate(Im{2,1}(:,:,z,c),reg_angle,'crop'); end
end
z_intensity{1}=squeeze(mean(mean(Y{1,3}(round(sz1/2)-100:round(sz1/2)+100,round(sz2/2)-100:round(sz2/2)+100,:))));
z_intensity{2}=squeeze(mean(mean(Y{2,3}(corr_coef(argm,1):corr_coef(argm,1)+200,corr_coef(argm,2):corr_coef(argm,2)+200,:))));
[cor cor_bin]=crosscorr(z_intensity{1},z_intensity{2},'numLags',size(Im{1,1},3)-1);
 [m argz]=sort(cor,'descend'); 
    g=1;
    while abs(cor_bin(argz(g)))>5
        g=g+1;
    end
    shift(1)=cor_bin(argz(g));
stack_shifted(1,:)=[1 size(Im{1,1},3)];    
stack_shifted(2,:)=[1-shift(1) size(Im{2,1},3)-shift(1)];
[stack_ini ini]=max(stack_shifted(:,1)); [stack_end endd]=min(stack_shifted(:,2));
for i=1:size(Y,1)
    for c=1:ch
    Y{i,c}=double(Y{i,c}(:,:,stack_ini-stack_shifted(i,1)+1:stack_end-stack_shifted(i,1)+1));
    end
    c_image_xy{i}=xyz_slice(Y{i,3});
end


 for i=1:2
ax{i} = axes('Units','pixels','Position',[50+450*(i-1) 50 400 400]);
imagesc(c_image_xy{i})
axis equal tight off 
hold all
title(D_name{i})
end

for dend=1:size(D_coor,2)
rng=round([min(D_coor{dend}); max(D_coor{dend})]);
  
[crop_im]=max(c_image_xy{1}(rng(1,2):rng(2,2),rng(1,1):rng(2,1),:),[],3);        
candidatePos{dend,1}=rng(1,1:2);
%draw_rectangle([rng(1,1:2) rng(2,1:2)-rng(1,1:2)+1],2,cmap(dend,:),ax{1});
edge=20;
[ssz1 ssz2]=size(max(c_image_xy{2},[],3));
clear search_im
search_im=zeros(ssz1+2*edge+1,ssz2+2*edge+1,1);
search_im(edge+1:edge+ssz1,edge+1:edge+ssz2)=max(c_image_xy{2},[],3);
[ match_data, ~ ] = matchPattern(search_im , crop_im,...
                                 0.1,2); 
pattRows = size(crop_im,1);  pattCols = size(crop_im,2);

%match_data(:,1) = match_data(:,1)+ceil(pattRows/2);  match_data(:,2) = match_data(:,2)+ceil(pattCols/2);  

[ss order]=sort(match_data(:,3),'ascend');
match_data=match_data(order,:);
match_data(:,1:2)=match_data(:,1:2)-edge;
candidatePos{dend,2}=[];
i=1;
while isempty(candidatePos{dend,2}) || i<=size(match_data,1) 
    
if sqrt(sum(([match_data(i,2) match_data(i,1)]-candidatePos{dend,1}).^2))<30
%     [m argm]=max(match_data(:,3));
%     match_data=match_data(argm,:);
    candidatePos{dend,2} =[match_data(i,2) match_data(i,1)];    
end
i=i+1;
end
if isempty(candidatePos{dend,2})
    candidatePos{dend,2}=[NaN NaN];
end
end
candidatePos=cell2mat(candidatePos);
NA=find(sum(isnan(candidatePos),2)==0);
for dend=NA'
plot(ax{1},D_coor{dend}(:,1),D_coor{dend}(:,2),'marker','.','color',cmap(dend,:),'linewidth',2,'markersize',18)
text(ax{1},D_coor{dend}(1,1),D_coor{dend}(1,2),num2str(dend),'color',cmap(dend,:))
s=[candidatePos(dend,3)-candidatePos(dend,1);
     candidatePos(dend,4)-candidatePos(dend,2)];
 plot(ax{2},D_coor{dend}(:,1)+s(1),D_coor{dend}(:,2)+s(2),'marker','.','color',cmap(dend,:),'linewidth',2,'markersize',18)
 text(ax{2},D_coor{dend}(1,1)+s(1),D_coor{dend}(1,2)+s(2),num2str(dend),'color',cmap(dend,:))
end

 end

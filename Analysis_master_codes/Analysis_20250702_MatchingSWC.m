f=88;
load(fullfile(fpath{f}, 'OP_Result.mat'), 'Result')%%
if ~isnan(StructureData{f})
folderPath = fileparts(StructureData{f});  % Change to your folder
else
    folderPath = fpath{f};  % Change to your folder
end
swcpath = dir(fullfile(folderPath, '*.swc'));
f
%%
ind=12; mag=1;
loadswc=fullfile(swcpath(ind).folder, swcpath(ind).name);
tree = load_tree(loadswc);
theta = 271;  % degrees
Rmat = [cosd(theta) -sind(theta) 0;
     sind(theta)  cosd(theta) 0;
     0            0           1];
Rmat=Rmat.*[mag mag 0;mag mag 0; 0 0 1];
[nROI nTime]=size(Result.normTraces);
sz=size(Result.ref_im); ftprnts=Result.ftprnt; cmap_ROIs=gen_colormap(Plasma,nROI);
if 1%isfield(Result,'tform')
RegT=affine2d(Rmat * Result.tform.T);
%RegT=affine2d(Rmat);
sx = sqrt(RegT.T(1,1)^2 + RegT.T(2,1)^2);  % scale along X
sy = sqrt(RegT.T(1,2)^2 + RegT.T(2,2)^2);  % scale along Y
str=[[transformPointsForward(RegT,[tree.X tree.Y])] tree.Z (tree.D+3)/mean([sx sy])];
else
    str=[tree.X tree.Y tree.Z tree.D+3];
end
somacoord=get_coord(Result.ftprnt(:,:,1));
offset=str(1,[1 2])-somacoord-[-2 3];
str(:,[1 2])=str(:,[1 2])-offset;

figure(19); clf; %show footprint
imshow2(Result.ref_im,[]); hold all;
scatter(str(:,1),str(:,2),str(:,4)+5,[0.6 0.6 0.6],'filled'); axis equal tight off; 
for r=1:nROI
DMDmask=ftprnts(:,:,r)>0;
str_ROI=str(find(str(:,1)>0.5 & str(:,1)<(sz(2)-0.5) & str(:,2)>0.5 & str(:,2)<(sz(1)-0.5)),:);
str_ROI=str_ROI(DMDmask(sub2ind(sz,round(str_ROI(:,2)),round(str_ROI(:,1)))),:);
scatter(str_ROI(:,1),str_ROI(:,2),str_ROI(:,4)/5,cmap_ROIs(r,:),'filled');
ROIpoly=bwboundaries(Result.ftprnt(:,:,r)>0);
plot(ROIpoly{1}(:,2),ROIpoly{1}(:,1),'color',cmap_ROIs(r,:),'LineWidth',1)
end
set(gca,'ydir','reverse'); hold all;
%bar100um=100/PCresult{f}.pixelsize;
%plot(prctile(str(:,1),95)-[0 bar100um],prctile(str(:,2),90)*[1 1],'color',[0 0 0],'LineWidth',3)
%%
Result.SWC=str;
Result.ParentMat=tree.dA;
dz=0.5;
Result.SWC=[Result.SWC(:,1:2) str(:,3)*dz Result.SWC(:,3)];
figure(20); clf;
scatter3(Result.SWC(:,1),Result.SWC(:,2),Result.SWC(:,3),Result.SWC(:,4)+5); axis equal; drawnow;

targetFile = strrep(fullfile(fpath{f},'OP_Result.mat'), '/Volumes/BHL18TB_D2/', '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/');
save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3')
save(targetFile,'Result','-v7.3')
disp(['Result is saved to ' targetFile]);
%%
% %%
% figure(22); clf;
% for f=foi
%     load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result')%%
%     nexttile([1 1]);
%     cmap_ROIs=gen_colormap(Plasma,size(Result.ftprnt,3));
%     ftprnt=Result.ftprnt(:,:,Result.dist_order);
%     sz=size(Result.ref_im); 
%     str=Result.SWC(:,[1:3]);
%     str(1,3)=50;
%     scatter(str(:,1),str(:,2),str(:,3)+15,[0.6 0.6 0.6],'filled'); axis equal tight off; hold all;
%     for r=1:size(Result.ftprnt,3)
% DMDmask=ftprnt(:,:,r)>0;
% str_ROI=str(find(str(:,1)>0.5 & str(:,1)<(sz(2)-0.5) & str(:,2)>0.5 & str(:,2)<(sz(1)-0.5)),:);
% str_ROI=str_ROI(DMDmask(sub2ind(sz,round(str_ROI(:,2)),round(str_ROI(:,1)))),:);
% scatter(str_ROI(:,1),str_ROI(:,2),str_ROI(:,3)+15,cmap_ROIs(r,:),'filled');
%     end
%     set(gca,'ydir','reverse'); hold all;
%     if f>=20
%         set(gca,'xdir','reverse');
%     end
% bar100um=100/PCresult{f}.pixelsize;
% plot(prctile(str(:,1),95)-[0 bar100um],prctile(str(:,2),90)*[1 1],'color',[0 0 0],'LineWidth',3);
% drawnow;
% end
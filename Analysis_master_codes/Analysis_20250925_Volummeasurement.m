clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Q192');

[~, ~, raw2] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:Z31');

fpath2read=[]; Str2read=[]; fpath2read_mouse=[];
ind2recon=find(cellfun(@length,raw(:,10))>20);
[~, unque_ind]=unique(cell2mat(raw(ind2recon,[2 5])),'rows');
fpath2read=raw(ind2recon(unque_ind),1);
fpath2read_mouse=raw(ind2recon(unque_ind),2);
Str2read=raw(ind2recon(unque_ind),10);
ind2recon=find(cellfun(@length,raw2(:,8))>20);
fpath2read=[fpath2read; raw2(ind2recon,1)];
fpath2read_mouse=[fpath2read_mouse; raw2(ind2recon,2)];
fpath2read_mouse=cell2mat(fpath2read_mouse);
Str2read=[Str2read; raw2(ind2recon,8)];
%%
%save(fullfile('/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/FigureS_VolumeReconstruct','Psf.mat'),"psfimg",'zprofile','roi','-v7.3')
load(fullfile('/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/FigureS2_VolumeReconstruct','Psf.mat'))
%load(fullfile('/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/FigureS_ScatteringMeasurement','STAmovs.mat'))
% psfimg=readtiff('/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/FigureS_VolumeReconstruct/Prism_20nmBead_05umStep_40um_10x.tif');
zcoord=[1:size(psfimg,3)]*0.5; %um
pfit=[];
for d=1:size(zprofile,2)
    [pfit(d,:) yhat R2]=fit_gauss1d(zcoord',zprofile(:,d));
end
sigma=mean(pfit(:,3),1)/0.42/0.42;
%%
fileList=[]; g=1; g2=1;
for f=1:length(fpath2read)
    f
    try
        load(fullfile(fpath2read{f},'OP_Result.mat'))
    catch
        load(fullfile(fpath2read{f},'PC_Result.mat'))
    end

    if isfield(Result,'SWC') 
        FtprntMask=max(Result.ftprnt,[],3)>0;
        [FtprntMask_indr FtprntMask_indc]=find(FtprntMask);
        if size(Result.SWC,2)==4
        fileList(g,:)=[f];
        Fraction_inplane=[];
        for z=1:max(Result.SWC(:,3)) %in um
            Fraction_inplane(z,:)=[length(find(Result.SWC(:,3)>z &  Result.SWC(:,3)<z+2*sigma))/size(Result.SWC,1) z+sigma];
        end

        [~, focusplane]=max(Fraction_inplane(:,1));
        InfocusIndex=[(Fraction_inplane(focusplane,2)-sigma)<Result.SWC(:,3) & (Fraction_inplane(focusplane,2)+sigma)>Result.SWC(:,3)...
            & ismember(round(Result.SWC(:,1)),FtprntMask_indc) & ismember(round(Result.SWC(:,2)),FtprntMask_indr)];
        InfocusVolume(g,:)=[sum(InfocusIndex)/size(Result.SWC,1)];
        g=g+1;
        end
        InFOVindex=[ismember(round(Result.SWC(:,1)),FtprntMask_indc) & ismember(round(Result.SWC(:,2)),FtprntMask_indr)];
        InFOVVolume(g2,:)=[sum(InFOVindex)/size(Result.SWC,1) fpath2read_mouse(f)];
        g2=g2+1;
    end
end

disp(['mean: ' num2str(mean(InFOVVolume(:,1)),3) ', std: ' num2str(std(InFOVVolume(:,1)),3)])

%%
repneuron=[2 3 4 12]; %  3 4 12
figure(28); clf; tiledlayout(2,2);
figure(29); clf; tiledlayout(4,1); pixsize=[0.936 1.17 1.17 0.936];
ax1=[]; g=1;
for f=fileList(repneuron)'
    try
        load(fullfile(fpath2read{f},'OP_Result.mat'))
    catch
        load(fullfile(fpath2read{f},'PC_Result.mat'))
    end

    Fraction_inplane=[];
    for z=1:max(Result.SWC(:,3))
        Fraction_inplane(z,:)=[length(find(Result.SWC(:,3)>z &  Result.SWC(:,3)<z+2*sigma))/size(Result.SWC,1) z+sigma];
    end

    [~, focusplane]=max(Fraction_inplane(:,1));
    FtprntMask=max(Result.ftprnt,[],3)>0;
    [FtprntMask_indr FtprntMask_indc]=find(FtprntMask);
    InfocusIndex=[(Fraction_inplane(focusplane,2)-sigma)<Result.SWC(:,3) & (Fraction_inplane(focusplane,2)+sigma)>Result.SWC(:,3)...
        & ismember(round(Result.SWC(:,1)),FtprntMask_indc) & ismember(round(Result.SWC(:,2)),FtprntMask_indr)];
    
    % StructureStack=(double(tiffreadVolume(Str2read{f})));
    % pixscale=[0.468 0.468 0.5];
    % x=[1:size(StructureStack,2)]*pixscale(1);
    % y=[1:size(StructureStack,1)]*pixscale(2);
    % z=[1:size(StructureStack,3)]*pixscale(3);

    % nexttile([1 1])
    % cax=[prctile(tovec(max(StructureStack,[],3)),15) prctile(tovec(max(StructureStack,[],3)),99.99)];
    % pcolor(x,y,max(StructureStack,[],3));
    % caxis(cax);
    % shading flat; axis equal tight off;
    % nexttile([1 1])
    % pcolor(y,z,max(permute(StructureStack,[1 3 2]),[],3)');
    % caxis(cax);
    % shading flat; axis equal tight off;
    % nexttile([1 1])
    % pcolor(x,z,max(permute(StructureStack,[2 3 1]),[],3)');
    % caxis(cax);
    % shading flat; axis equal tight off;
    figure(28);
    nexttile([1 1])
    sizevec=ones(1,size(Result.SWC,1))*2; sizevec(1)=30;
    scatter3((max(Result.SWC(:,1))-Result.SWC(:,1))*pixsize(g),Result.SWC(:,2)*pixsize(g),Result.SWC(:,3),sizevec,[0 0 0],'filled'); axis equal; hold all;
    %scatter3(Result.SWC(InfocusIndex,1)*1.17,Result.SWC(InfocusIndex,2)*1.17,Result.SWC(InfocusIndex,3)*0.5,sizevec(InfocusIndex),[1 0 0],'filled'); axis equal;
    %xlabel('Distance (\mum)'); ylabel('Distance (\mum)'); zlabel('Distance (\mum)');
    zl = zlim;
    zticks(round(linspace(zl(1), zl(2), 3)));
    set(gca,'YDir','reverse')
    view([210-180 30])
    figure(29);
    ax1=[nexttile([1 1])];
    %scatter(max(Result.SWC(:,1)*pixsize(g))-Result.SWC(:,1)*pixsize(g),Result.SWC(:,3),sizevec,[0 0 0],'filled'); hold all; axis equal tight; grid on;
    scatter(Result.SWC(:,1)*pixsize(g),Result.SWC(:,3),sizevec,[0 0 0],'filled'); hold all; axis equal tight; grid on;
    fill([0 max(Result.SWC(:,1)*pixsize(g)) max(Result.SWC(:,1)*pixsize(g)) 0],Fraction_inplane(focusplane,2)+[-sigma -sigma sigma sigma],[1 0 0],'EdgeColor','none','FaceAlpha',0.3)
    xlabel('Distance (\mum)'); ylabel('Distance (\mum)');
    %scatter(Result.SWC(InfocusIndex,1)*1.17,Result.SWC(InfocusIndex,3)*0.5,sizevec(InfocusIndex),[1 0 0],'filled'); axis equal;
    %view([210 0])
    g=g+1;
end
set_fontsize(13);
figure(28);
set_fontsize(13);
%linkObj = linkprop([ax1],{'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});



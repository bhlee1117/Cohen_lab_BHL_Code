% Analysis on AAV expression sample and plot, in house YQ201
% 

clear
[fpath] = uigetfile_n_dir;
%%
for i=[1:length(fpath)]
mov=vm(fpath{i});
load(fullfile(fpath{i},'settings.mat'));
mov_test=mov(:,:,150:250);
try mov_test = single(mov_test)./single(max(mov_test.data(:)));
catch disp('change to vm')

mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,xyField]=optical_flow_motion_correction_LBH((mov),mov_ref,'optic_flow');
clear mov
mov_mc=vm(mov_mc);
mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
mcTrace = squeeze(mean(xyField,[1 2]));
save([fpath{i} '/mcTrace.mat'],'mcTrace')
clear mov_mc
delete([fpath{i} '/Sq_camera.bin'])
end
%load('Result_221107.mat')

%%
R=[14];
for i=1:length(fpath)
    Sz = importdata([fpath{i} '/experimental_parameters.txt']);
    sz1=Sz.data(1); sz2=Sz.data(2);
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz2,sz1));
    im_mean{i}=mean(mov_mc,3);
    im_G=imgaussfilt(mean(mov_mc,3),2);
    [centers radii]=Cell_segment_circle_10x(im_G);
    centers=cell_detection_manual(mean(mov_mc,3),centers);
    mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
    for j=1:length(R)
        c_ftprnt{i,j}=mask_footprint(centers,mov_res,[],R(j));
        traces_bin{i,j}=-(tovec(mov_mc)'*tovec(c_ftprnt{i,j}>0))';
        traces{i,j}=-(tovec(mov_mc)'*tovec(c_ftprnt{i,j}))';

        t=traces{i,j}-movmean(traces{i,j},200,2);
        t2=traces{i,j}-movmean(traces{i,j},30,2);
        clear noise
        for n=1:size(traces{i,j},1)
            noise(n)=get_threshold(t2(n,:),1); 
        end
        spike{i,j}=find_spike_bh(t./noise',4);
        SNR_tr{i,j}=t./noise';
    end
end

%%
i=3;
show_footprnt(c_ftprnt{i,1},im_mean{i})
show_traces(SNR_tr{i,1},10)
%%
g=1; L=[]; foi=[1 2 3 4]; n=[35 35 36 37; 37 37 33 34; 33 33 29 30];
tiledlayout(size(c_ftprnt,2)+1,length(foi));
for i=foi
    nexttile([1 1])
    imshow2(im_mean{i},[])
    %title(['Target illumination sq size: ' num2str(g*2)])

    g=g+1;
end
g=1; gg=1;
for j=1:size(c_ftprnt,2)
    gg=1;
    for i=foi

        nexttile([1 1])
        imshow2(max(c_ftprnt{i,g}(:,:,n(:,gg)),[],3),[])
        %imshow2(max(c_ftprnt{i,g},[],3),[])
        %title(['Footprint radius: ' num2str(g*3)])
        gg=gg+1;
    end
    L=[L {['R = ' num2str(g*3)]}];
    g=g+1;
end
%%
clear SNR_tr spike
g=1; foi=[1 2]; cmap=distinguishable_colors(6);
figure; tiledlayout(2,3);
ax=[]; n=[35 35 36 37; 37 37 33 34; 33 33 29 30];
for i=foi %file to open
    li=[];
    [ident_list{g} match{g}]=match_cell(c_ftprnt(i,:),5);
    ax=[ax nexttile([1 1])];
    for j=1:size(c_ftprnt,2) % refer footprint
        %for c=size(ident_list{g},1) % cell
        for c=1:size(ident_list{g},1) % cell
            t=traces{i,j}(ident_list{g}(c,j),:)-movmean(traces{i,j}(ident_list{g}(c,j),:),200,2);
            t2=traces{i,j}(ident_list{g}(c,j),:)-movmean(traces{i,j}(ident_list{g}(c,j),:),30,2);
            noise=get_threshold(t2,1);
            spike{i,j}(ident_list{g}(c,j),:)=find_spike_bh(t,noise*4);
            SNR_tr{i,j}(ident_list{g}(c,j),:)=t/noise;
            li=[li plot(SNR_tr{i,j}(ident_list{g}(c,j),:),'color',cmap(j,:))];
            hold all
            plot(find(spike{i,j}(ident_list{g}(c,j),:)),t(find(spike{i,j}(ident_list{g}(c,j),:)))/noise,'r.')
        end
    end
    linkaxes(ax,'xy')
    g=g+1;

end
%lgd=legend(li,L); legend('boxoff');
%lgd.Layout.Tile = 'east';

%% Plot SNR
clear dat
g=1;
tiledlayout(1,length(foi));
ax=[];
for i=foi
    [ident_list{g} match{g}]=match_cell(c_ftprnt(i,:),5);
    ax=[ax nexttile([1 1])];
    s=sum(cell2mat(spike(i,:)'),1);
    for j=1:size(c_ftprnt,2)
        dat{g}(:,j)=SNR_tr{i,j}(ident_list{g}(end,j),find(s>0));
    end
    %plot([1:5]*3+randn(size(dat{g},1),5)*0.1,dat{g},'.','color',[0.6 0.6 0.6])
    plot(repmat([1:5]'*3,1,size(dat{g},1)),dat{g}')
    hold all
    errorbar([1:5]*3,mean(dat{g},1),std(dat{g},0,1),'k')
    ylabel('SNR'); xlabel('Footprint radius');
    xlim([1.5 16.5])
    g=g+1;
end
linkaxes(ax,'y')
%% High-pass filter
mov_mc=double(mov);
im_G=imgaussfilt(mean(mov_mc,3),2);
mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[8 8]),2);
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[16 16]),2);
[centers radii]=Cell_segment_circle_10x(im_G);
centers=cell_detection_manual(mean(mov_mc,3),centers);
%%
mcTrace = squeeze(mean(xyField,[1 2]));
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));
c_ftprnt=mask_footprint(centers,mov_res,[],10);
show_footprnt(c_ftprnt,mov_mc)
coord=get_coord(c_ftprnt);

traces=-(tovec(mov_mc)'*tovec(c_ftprnt>0))';
traces_hi=squeeze(traces) - movmean(squeeze(traces),400,2);
%traces_hi=squeeze(traces) - mean(squeeze(traces),2);
traces_hi=traces_hi./range(traces_hi,2);
%traces_hi=-traces_hi./mean(traces_hi,2);
%%
figure; scale=0.6;
tiledlayout(10,4)
ax1 = nexttile([2 2]);
colr = max(colormap(jet(size(c_ftprnt,3))),0);
imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

ax4 = nexttile([2 2]);
imshow2(mean(mov_mc,3),[]); colormap('gray')

linkaxes([ax1 ax4],'xy')

ax2 = nexttile([7 4]);
lines=plot(traces_hi'+[1:size(traces_hi,1)]*scale);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(size(c_ftprnt,3)),2))
axis tight
set(gca,'ytick',[1:size(traces_hi,1)]*scale,'yticklabel',[1:size(traces_hi,1)])
ax3 = nexttile([1 4]);
plot(Blue)

linkaxes([ax2 ax3],'x')
saveas(gcf,[fpath 'voltage_trace_plot'])
save([fpath 'result.mat'],'c_ftprnt','traces')

function dendrite_image_show(Y3,window_sz,weight)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 window_sz*6 1000]);
avail=find(sum(cell2mat(cellfun(@isempty,Y3,'UniformOutput',false)))==0);
cmap=[1 0 0; 0 1 1;1 0 1;1 1 0];
for dend=avail
    clf('reset')
    
    for d=1:2
        cc=find(max(max(Y3{d,dend},[],3),[],1)>0);
        rr=find(max(max(Y3{d,dend},[],3),[],2)>0);
        minC=min(cc); maxC=max(cc); minR=min(rr); maxR=max(rr);
        im_show=zeros(maxR-minR+1,maxC-minC+1,3);
        for ch=1:4
            each_chan_im_show{ch}=zeros(maxR-minR+1,maxC-minC+1,3);
            for c=1:3
                t=sort(unique(Y3{d,dend}(:,:,ch)),'ascend');
                im_show(:,:,c)=im_show(:,:,c)+(Y3{d,dend}(minR:maxR,minC:maxC,ch)-t(round(length(t)*weight(ch)),1))/t(end)*cmap(ch,c);
                each_chan_im_show{ch}(:,:,c)=each_chan_im_show{ch}(:,:,c)+(Y3{d,dend}(minR:maxR,minC:maxC,ch)-t(round(length(t)*weight(ch)),1))/t(end)*cmap(ch,c);
            end
            ax{d,ch} = axes('Units','pixels','Position',[50+window_sz*(ch) 50+450*(d-1) 400 400]);
            [r c]=find(max(im_show,[],3)>0);
            if size(each_chan_im_show{ch},1)<size(each_chan_im_show{ch},2)
            imagesc(imrotate(each_chan_im_show{ch},90,'loose'))
            else
            imagesc(each_chan_im_show{ch})    
            end
            axis equal tight off
            hold all
            title(['Channel', num2str(ch)])
        end
        ax{d,ch} = axes('Units','pixels','Position',[50 50+450*(d-1) 400 400]);
        [r c]=find(max(im_show,[],3)>0);
        if size(im_show,1)<size(im_show,2)
            imagesc(imrotate(im_show,90,'loose'))
            else
        imagesc(im_show)
            end
        
        axis equal tight off
        hold all
        title(['Dendrite #',num2str(dend),' on day ',num2str(d)])
        
    end
    waitforbuttonpress
end
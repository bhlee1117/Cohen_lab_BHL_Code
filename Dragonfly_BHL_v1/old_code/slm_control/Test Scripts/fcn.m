            imgbg = prctile(newmov(:),95);
            close all
            moms = imgmoments(max(newmov-imgbg,0)); 
            t = 0:1:360;
            for itSortStep = 1:nSteps
                figure windowstyle docked
                imshow(newmov(:,:,itSortStep),[0 5e3])
                hold on
%                 ar = -ihs:ihs;
%                 [ay, ax] = ndgrid(ar,ar); % grid
%                 fitmat = [ones(size(ay(:))) ax(:) ay(:) ax(:).^2 ay(:).^2 ax(:).*ay(:)];
%                 imgdata = reshape(newmov(:,:,itSortStep),[],1);
%                 x = (imgdata.^2.*fitmat)\(imgdata.^2.*log(imgdata));
%                 
% 
%                 cx = mean(mean(newmov(:,:,itSortStep).*(1:121),1),2)./mean(mean(newmov(:,:,itSortStep),1),2);
%                 cy = mean(mean(newmov(:,:,itSortStep).*(1:121)',1),2)./mean(mean(newmov(:,:,itSortStep),1),2);
%                 sx = sqrt(mean(mean(newmov(:,:,itSortStep).*((1:121)-cx).^2,1),2)./mean(mean(newmov(:,:,itSortStep),1),2));
%                 sy = sqrt(mean(mean(newmov(:,:,itSortStep).*((1:121)'-cy).^2,1),2)./mean(mean(newmov(:,:,itSortStep),1),2));
                plot(ihs+1,ihs+1,'ob')
                plot(ihs+1+moms(4,itSortStep),ihs+1+moms(5,itSortStep),'or')
                plot(...
                    ihs+1+moms(4,itSortStep)+moms(6,itSortStep)*cosd(t),...
                    ihs+1+moms(5,itSortStep)+moms(7,itSortStep)*sind(t))
%                 plot(cx,cy,'o')
%                 plot(...
%                     cx+sx*cosd(t),...
%                     cy+sy*sind(t))
            end
        

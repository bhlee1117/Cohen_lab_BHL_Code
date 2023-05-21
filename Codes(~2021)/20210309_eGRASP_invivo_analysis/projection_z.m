function Y3=projection_z(Y2,coor,Pos,pix_xyz,range)
avail=find(sum(isnan(Pos),2)==0);
xz_rat=pix_xyz(1,2)/pix_xyz(1,1);
cmap=distinguishable_colors(size(Pos,1));

g=1; cmap=[1 0 0;0 1 1;1 0 1;1 1 0];

for dend=avail'
coor{dend}=coor{dend}(find(sum(isnan(coor{dend}),2)==0),:);
if size(coor{dend},1)>2
    clear V V2 V3 
    vs=(inv(coor{dend}([1 round((size(coor{dend},1)+1)/2) size(coor{dend},1)],:))*[1;1;1])';
    vs=vs/sqrt(sum(vs.^2));
    
    s=[Pos(dend,3)-Pos(dend,1)  Pos(dend,4)-Pos(dend,2) 0];
    if 0%max(coor{dend}(:,3))-min(coor{dend}(:,3))>15
%         z_rng=([min(coor{dend}(:,3)):(max(coor{dend}(:,3))-min(coor{dend}(:,3)))/100:max(coor{dend}(:,3))]);
%         int_dend{1}=[interp1(coor{dend}(:,3),coor{dend}(:,1),z_rng,'linear')' interp1(coor{dend}(:,3),coor{dend}(:,2),z_rng,'linear')' z_rng'];
    else
        try %  max(coor{dend}(:,2))-min(coor{dend}(:,2))   > max(coor{dend}(:,1))-min(coor{dend}(:,1))
            y_rng=([min(coor{dend}(:,2)):(max(coor{dend}(:,2))-min(coor{dend}(:,2)))/100:max(coor{dend}(:,2))]);
            int_dend{1}=[interp1(coor{dend}(:,2),coor{dend}(:,1),y_rng,'linear')' y_rng' interp1(coor{dend}(:,2),coor{dend}(:,3),y_rng,'linear')'];
        catch
             [au,ia] = unique(round(coor{dend},1),'stable');
                Same = ones(size(coor{dend}));
                Same(ia) = 0;
                coor{dend}=coor{dend}(find(sum(Same,2)==0),:);
            x_rng=([min(coor{dend}(:,1)):(max(coor{dend}(:,1))-min(coor{dend}(:,1)))/100:max(coor{dend}(:,1))]);
            int_dend{1}=[x_rng' interp1(coor{dend}(:,1),coor{dend}(:,2),x_rng,'linear')' interp1(coor{dend}(:,1),coor{dend}(:,3),x_rng,'linear')'];
        end
    end
    int_dend{2}=int_dend{1}+s;
    for d=1:2
        [sz1 sz2 sz3]=size(Y2{d,3});
        [xx yy zz]=meshgrid([1:sz2],[1:sz1],[1:sz3]);
        V=zeros(sz1,sz2,sz3);
        for z=1:size(int_dend{d},1)-1
            V(sqrt((xx-int_dend{d}(z,1)).^2+(yy-int_dend{d}(z,2)).^2++xz_rat*(zz-int_dend{d}(z,3)).^2)<range)=1;
        end
        
        % V2=imgaussfilt(V,4);
        % p = patch(isosurface(xx,yy,zz,V2,0.8));
        % isonormals(xx,yy,zz,V2,p)
        % p.FaceColor = [1 1 0];
        % p.EdgeColor = 'none';
        % hold all
        
        for ch=1:size(Y2,2)
            mask_im=Y2{d,ch};
            mask_im(V==0)=0;
            Y3{d,dend}(:,:,ch)=max(imrotate3(mask_im,rad2deg(acos(sum(vs.*[0 0 1],2)/sqrt(sum((vs).^2)))),[vs(2) vs(1) 0]/sqrt(sum((vs).^2))),[],3);
            %V3=imrotate3(V2,rad2deg(acos(sum(vs.*[0 0 1],2)/sqrt(sum((vs).^2)))),[vs(2) vs(1) 0]/sqrt(sum((vs).^2)));
            
        end
        
    end
    % daspect([1 1 1])
    % view(3);
    % axis tight
    % camlight
    % lighting gouraud
    
    % [sz1 sz2 sz3]=size(V3);
    % [xx yy zz]=meshgrid([1:sz2],[1:sz1],[1:sz3]);
    % p = patch(isosurface(xx,yy,zz,V3,0.8));
    % isonormals(xx,yy,zz,V3,p)
    % p.FaceColor = [1 0 0];
    % p.EdgeColor = 'none';
    % hold all
    % daspect([1 1 1])
    % view(3);
    % axis tight
    % camlight
    % lighting gouraud
end
end
%%


end
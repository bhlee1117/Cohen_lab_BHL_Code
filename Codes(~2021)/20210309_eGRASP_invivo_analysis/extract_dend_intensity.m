function Dend_int=extract_dend_intensity(Y2,coor,Pos,pix_xyz,range)
avail=find(sum(isnan(Pos),2)==0);
xz_rat=pix_xyz(1,2)/pix_xyz(1,1);
cmap=distinguishable_colors(size(Pos,1));
for dend=avail'
    dend
    s=[Pos(dend,3)-Pos(dend,1)  Pos(dend,4)-Pos(dend,2) 0];
    if 0% max(coor{dend}(:,3))-min(coor{dend}(:,3))>20
        % z_rng=([min(coor{dend}(:,3)):(max(coor{dend}(:,3))-min(coor{dend}(:,3)))/100:max(coor{dend}(:,3))]);
        % int_dend{1}=[interp1(coor{dend}(:,3),coor{dend}(:,1),z_rng,'linear')' interp1(coor{dend}(:,3),coor{dend}(:,2),z_rng,'linear')' z_rng'];
    else
        try  %max(coor{dend}(:,2))-min(coor{dend}(:,2))   > max(coor{dend}(:,1))-min(coor{dend}(:,1))
            y_rng=([min(coor{dend}(:,2)):(max(coor{dend}(:,2))-min(coor{dend}(:,2)))/100:max(coor{dend}(:,2))]);
            int_dend{1}=[interp1(coor{dend}(:,2),coor{dend}(:,1),y_rng,'linear')' y_rng' interp1(coor{dend}(:,2),coor{dend}(:,3),y_rng,'linear')'];
        catch
            try
                x_rng=([min(coor{dend}(:,1)):(max(coor{dend}(:,1))-min(coor{dend}(:,1)))/100:max(coor{dend}(:,1))]);
                int_dend{1}=[x_rng' interp1(coor{dend}(:,1),coor{dend}(:,2),x_rng,'linear')' interp1(coor{dend}(:,1),coor{dend}(:,3),x_rng,'linear')'];
            catch
                [au,ia] = unique(round(coor{dend},1),'stable');
                Same = ones(size(coor{dend}));
                Same(ia) = 0;
                coor{dend}=coor{dend}(find(sum(Same,2)==0),:);
                x_rng=([min(coor{dend}(:,1)):(max(coor{dend}(:,1))-min(coor{dend}(:,1)))/100:max(coor{dend}(:,1))]);
                int_dend{1}=[x_rng' interp1(coor{dend}(:,1),coor{dend}(:,2),x_rng,'linear')' interp1(coor{dend}(:,1),coor{dend}(:,3),x_rng,'linear')'];
            end
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
        % % p.FaceColor = cmap(dend,:);
        % p.EdgeColor = 'none';
        % hold all
        Dend_int{d,dend}=Y2{d,1}(V==1);
    end
end
% daspect([1 1 1])
% view(3);
% axis tight
% camlight
% lighting gouraud
end
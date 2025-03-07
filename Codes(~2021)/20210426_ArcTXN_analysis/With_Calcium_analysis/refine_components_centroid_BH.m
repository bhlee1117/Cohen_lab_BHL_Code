function [A,C,newcenters,CC] = refine_components_centroid_BH(Y,A,C,centers,img,sx,options,rm_centers,ad_centers)
% function that allows to add or remove components. 
% Usage: 
% after runnning use left click to remove and right click to add components
% at specific locations. Hit enter when done
% Inputs
% ------
% Y:       data loaded in memory
% A:       current set of spatial components (excluding background)
% C:       current set of temporal components (excluding background)
% centers: coordinates of the centers of the components
% img:     image used to identify components (for instance correlation image)
% sx:      half-width of search window
% options: options structure
% 
% Returns:
% --------
% A:          updated set of spatial components
% C:          updated set of temporal components
% newcenters: updated component centers
% CC:         contour plots of components
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 700 900]);
defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = defoptions; end
if ~isfield(options,'d1') || isempty(options.d1); options.d1 = input('What is the total number of rows? \n'); end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); options.d2 = input('What is the total number of columns? \n'); end       % # of columns
if ~isfield(options,'cont_threshold') || isempty(options.cont_threshold); cont_threshold = defoptions.cont_threshold; else cont_threshold = options.cont_threshold; end          % # of rows
if nargin < 6 || isempty(sx)
    sx = 5;
end
if nargin < 5 || isempty(img)
    img = std(Y,[],3);
end
if nargin < 4 || isempty(centers)
    centers = com(A,options.d1,options.d2);
end
    

min_distance_point_selection=2;
x=1;
[K,T] = size(C);
newcenters=[centers];
CC=[];
for i = 1:size(A,2)
    a_srt = sort(A(:,i),'descend');
    ff = find(cumsum(a_srt.^2) >= cont_threshold*sum(a_srt.^2),1,'first');
%     CC{i} = contour(reshape(A(:,i),options.d1,options.d2),[0,0]+a_srt(ff),'Linecolor',[1,0,1]/2);
%     CC{i}(CC{i}<1) = NaN;
%     CC{i}(:,CC{i}(1,:)>options.d2) = NaN;
%     CC{i}(:,CC{i}(2,:)>options.d1) = NaN;
    hold on;
end
    
% remove
for i=1:size(rm_centers,2)
            [y]=centers(rm_centers(1,i),1); [x]=centers(rm_centers(1,i),2);
            [m,id]=min(sum(bsxfun(@minus,newcenters,[y,x]).^2,2));            
            ident_point=[newcenters(id,2),newcenters(id,1)];
            if m<=min_distance_point_selection
                disp(['Removing point:' num2str(ident_point)]) 
                newcenters(id,:)=[]; 
                A(:,id) = [];
                C(id,:) = [];
%                 CC(id) = [];
                K = K - 1;
            end
end
               
            
% add            
for i=1:size(ad_centers,1)
[y]=ad_centers(i,1); [x]=ad_centers(i,2);
        pixel=round([x y]);
            newcenters=[newcenters; fliplr(pixel)];
            int_x = round(newcenters(end,1)) + (-sx:sx);
            if int_x(1)<1
                int_x = int_x + 1 - int_x(1);
            end
            if int_x(end)>options.d1
                int_x = int_x - (int_x(end)-options.d1);
            end
            int_y = round(newcenters(end,2)) + (-sx:sx);
            if int_y(1)<1
                int_y = int_y + 1 - int_y(1);
            end
            if int_y(end)>options.d2
                int_y = int_y - (int_y(end)-options.d2);
            end
            [INT_x,INT_y] = meshgrid(int_x,int_y);
            coor = sub2ind([options.d1,options.d2],INT_x(:),INT_y(:));
            Ypatch = reshape(Y(int_x,int_y,:),(2*sx+1)^2,T);                        
            Y_res = Ypatch - A(coor,:)*C;
            Y_res = bsxfun(@minus, Y_res, median(Y_res,2));
            [atemp, ctemp, ~, ~, newcenter, ~] = greedyROI(reshape(Y_res,2*sx+1,2*sx+1,T), 1, options);
            %[atemp, ctemp] = initialize_components(reshape(Y_res,2*sx+1,2*sx+1,T), 1,sx,options);  % initialize
            % find contour
            a_srt = sort(atemp,'descend');
            ff = find(cumsum(a_srt.^2) >= cont_threshold*sum(a_srt.^2),1,'first');
            K = K + 1;
            A(coor,K) = atemp/norm(atemp);
            C(K,:) = ctemp*norm(atemp);
            new_center = com(A(:,end),options.d1,options.d2);
            newcenters(end,:) = new_center;
%             scatter(new_center(2),new_center(1),'mo'); hold on; 
%             CC{K} = contour(reshape(A(:,end),options.d1,options.d2),[0,0]+a_srt(ff),'Linecolor',[1,0,1]/2);
%             CC{K}(CC{K}<1) = NaN;
%             CC{K}(:,CC{K}(1,:)>options.d2) = NaN;
%             CC{K}(:,CC{K}(2,:)>options.d1) = NaN;
end
end
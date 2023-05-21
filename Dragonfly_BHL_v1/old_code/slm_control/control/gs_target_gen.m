function Target = gs_target_gen(roi,type)
switch type
    case 'spt'
        spt_sz = 1;
        gauss2 = @(x,y,x0,y0,sigma) exp(-1/sigma^2/2*((x-x0).^2+(y-y0).^2))/2/pi/sigma^2;

        [Y,X] = ndgrid(1:1152,1:1920);
        Target = zeros(1152,1920);
        % spot = zeros(1152,1920);
        % spot(round(576-spt_sz/2/n):round(576+spt_sz/2/n),round(960-spt_sz/2):round(960+spt_sz/2))=1;
        % theta = linspace(2*pi/10,2*pi,10);
        if iscell(roi)
            for i=1:length(roi)
                x_spt(i) = roi{i}(1);
                y_spt(i) = roi{i}(2);
            end
        else
            x_spt = roi(:,1);
            y_spt = roi(:,2);
        end
        for i = 1:length(x_spt)

        Target = Target + gauss2(X,Y,x_spt(i),y_spt(i),spt_sz);
        end
    case 'roi'
        [X,Y]=ndgrid(1:1152,1:1920);

        in = cellfun(@(pixel_pos) inpolygon(X,Y,pixel_pos(:,2),pixel_pos(:,1)),roi,'uniformoutput',false);
        Target = false(1152,1920);
        for i=1:length(in), Target = Target|in{i}; end
    case 'roi2spt'
        [X,Y]=ndgrid(1:1152,1:1920);

        in = cellfun(@(pixel_pos) inpolygon(X,Y,pixel_pos(:,2),pixel_pos(:,1)),roi,'uniformoutput',false);
        x_spt = [];
        y_spt = [];    
        area = 10*20;
        for i=1:length(in)
            in_i = find(in{i});
            [x_spt_i,y_spt_i] = ind2sub(size(in{i}),in_i(1:area:end));
            x_spt = [x_spt; x_spt_i];
            y_spt = [y_spt; y_spt_i];
        end
        Target = [x_spt y_spt];
%         spt_sz = 1;
%         gauss2 = @(x,y,x0,y0,sigma) exp(-1/sigma^2/2*((x-x0).^2+(y-y0).^2))/2/pi/sigma^2;
% 
%         [X,Y] = ndgrid(1:1152,1:1920);
%         Target = zeros(1152,1920);
% 
%         for i = 1:length(x_spt)
%             Target = Target + gauss2(X,Y,x_spt(i),y_spt(i),spt_sz);
%         end
end
end
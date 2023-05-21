%% square 1k
% 	 Target = ones(1000,1000);
%     m=2;
%     n=8;
%     Target = padarray(imresize(Target,floor([1152/n 1920]/m)),[ceil(1152*(1-1/m*1/n)/2) ceil(1920*(1-1/m)/2)],'both');
%     Target = imresize(Target,[1152 1920]);
%% 
% 
% %m=2;
% n=2;
% 
% spt_sz = 1;
% gauss2 = @(x,y,x0,y0,sigma) exp(-1/sigma^2/2*((x-x0).^2+(y-y0).^2))/2/pi/sigma^2;
% 
% [X,Y] = ndgrid(-1152/2:1152/2-1,-1920/2:1920/2-1);
% Target = zeros(1152,1920);
% % spot = zeros(1152,1920);
% % spot(round(576-spt_sz/2/n):round(576+spt_sz/2/n),round(960-spt_sz/2):round(960+spt_sz/2))=1;
% theta = linspace(2*pi/5,2*pi,5);
% % for i = 1:length(theta)
% 
% % % Target = Target + gauss2(X,Y,cos(theta(i))*100,sin(theta(i))*100,spt_sz);
% % end
% 
% Target = gauss2(X,Y,0,0,spt_sz);
% % Target(1152-(800:850),700:1400)=1;
% % Target = Target + spot;
% % 
% % figure(33);clf;imshow(Target,[])
%%
% Target = register_pattern_test(1152,1920,[700 500]);
% figure(10);clf;imshow(Target,[])




%% clicky
% test_im = zeros(1152,1920);
% [ROIs,~]=clicky_faster(test_im);
% [X,Y]=ndgrid(1:1152,1:1920);
% 
% in = cellfun(@(pixel_pos) inpolygon(X,Y,pixel_pos(:,2),pixel_pos(:,1)),ROIs,'uniformoutput',false);
% Target = false(1152,1920);
% for i=1:length(in), Target = Target|in{i}; end

%% click spots
% test_im = zeros(1152,1920);
% % [ROIs,~]=clicky_faster(test_im);
% 
% pts = click_im(test_im,'pts');
%  
% x_spt = pts(:,1);
% y_spt = pts(:,2);
% 
% spt_sz = 1;
% gauss2 = @(x,y,x0,y0,sigma) exp(-1/sigma^2/2*((x-x0).^2+(y-y0).^2))/2/pi/sigma^2;
% 
% [Y,X] = ndgrid(1:1152,1:1920);
% Target = zeros(1152,1920);
% % spot = zeros(1152,1920);
% % spot(round(576-spt_sz/2/n):round(576+spt_sz/2/n),round(960-spt_sz/2):round(960+spt_sz/2))=1;
% % theta = linspace(2*pi/10,2*pi,10);
% for i = 1:length(x_spt)
% 
% Target = Target + gauss2(X,Y,x_spt(i),y_spt(i),spt_sz);
% end
%% clicky into discrete spots
test_im = zeros(1152,1920);
% [ROIs,~]=clicky_faster(test_im);
ROIs = click_im(test_im,'rect');

[X,Y]=ndgrid(1:1152,1:1920);

in = cellfun(@(pixel_pos) inpolygon(X,Y,pixel_pos(:,2),pixel_pos(:,1)),ROIs,'uniformoutput',false);
x_spt = [];
y_spt = [];    
area = 25*25;
for i=1:length(in)
    in_i = find(in{i});
    [x_spt_i,y_spt_i] = ind2sub(size(in{i}),in_i(1:area:end));
    x_spt = [x_spt; x_spt_i];
    y_spt = [y_spt; y_spt_i];
end

spt_sz = 1;
gauss2 = @(x,y,x0,y0,sigma) exp(-1/sigma^2/2*((x-x0).^2+(y-y0).^2))/2/pi/sigma^2;

[X,Y] = ndgrid(1:1152,1:1920);
Target = zeros(1152,1920);
% spot = zeros(1152,1920);
% spot(round(576-spt_sz/2/n):round(576+spt_sz/2/n),round(960-spt_sz/2):round(960+spt_sz/2))=1;
% theta = linspace(2*pi/10,2*pi,10);
for i = 1:length(x_spt)

Target = Target + gauss2(X,Y,x_spt(i),y_spt(i),spt_sz);
end

% Target = gauss2(X,Y,0,0,spt_sz)+gauss2(X,Y,20,20,spt_sz);
% Target(1152-(800:850),700:1400)=1;
% Target = Target + spot;

% figure(33);clf;imshow(Target,[])
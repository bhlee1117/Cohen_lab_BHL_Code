function [x_spt,y_spt] = get_pts()

test_im = zeros(1152,1920);
% [ROIs,~]=clicky_faster(test_im);

pts = click_im(test_im,'pts');
 
x_spt = pts(:,1);
y_spt = pts(:,2);

end
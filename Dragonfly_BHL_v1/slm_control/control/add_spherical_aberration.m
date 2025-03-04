


cm = hSLMApp.current_mask;

[xx, yy] = ndgrid(1:size(cm,1),1:size(cm,2));
xx = (xx-1920/2)/1920;
yy = (yy-1152/2)/1152;
rr = sqrt(xx.^2+yy.^2);
th = atan2(yy,xx);

sa = sqrt(5)*(6*rr.^4-6*rr.^2+1); % Primary spherical

cm_sa = mod(cm+sa*5,2*pi);

hSLM.project(cm_sa)
function phase = wgs_spots(spt_corrd, iter, isGPU)

f = 1e-1; % meter effective focal length (guess)
d = 9.2e-6; % SLM pixel size
lambda = .594e-6;

H = 1152; L = 1920; % SLM dimension

slm_mask = false(1152, 1920);
slm_mask(1152/4+1:1152*3/4,1920/4+1:1920*3/4)=true;

% SLM face light intensity field
x = linspace(-15,15,1920);
y = linspace(-10,10,1152);
[X,Y] = meshgrid(x,y);
x0 = 0;     		% center
y0 = 0;     		% center
sigma = 2; 			% beam waist
A = 1;      		% peak of the beam 
res = ((X-x0).^2 + (Y-y0).^2)./(2*sigma^2);
input_intensity = A  * exp(-res).*slm_mask;
input_intensity = input_intensity/sqrt(sum(input_intensity.^2,[1 2]));

x = reshape(spt_corrd(:,1),1,1,[])*d;
y = reshape(spt_corrd(:,2),1,1,[])*d;
z = reshape(spt_corrd(:,3),1,1,[])*d;

[xj,yj] = meshgrid(((1:L)-L/2)*d,((1:H)-H/2)*d);

wt = ones(1,1,length(x));

% [xx, yy] = ndgrid(1:size(xj,1),1:size(xj,2));
% xx = (xx-1920/2)/1920;
% yy = (yy-1152/2)/1152;
% rr = sqrt(xx.^2+yy.^2);
% th = atan2(yy,xx);
% 
% sa = sqrt(5)*(6*rr.^4-6*rr.^2+1); % Primary spherical

% sptPh = pi*z/lambda/f^2.*(xj.^2+yj.^2)+2*pi/f/lambda*(xj.*x+yj.*y)+sa*10;
sptPh = pi*z/lambda/f^2.*(xj.^2+yj.^2)+2*pi/f/lambda*(xj.*x+yj.*y);
sptPhOffset = rand(size(xj))*2*pi;
slmField = sum(exp(1i*(-sptPhOffset-sptPh)),3);

if isGPU
    wt = gpuArray(wt);
    sptPh = gpuArray(sptPh);
    sptPhOffset = gpuArray(sptPhOffset);
    slmField = gpuArray(slmField);
%     slmPh = gpuArray(slmPh);
end

for i=1:iter
    slmPh = angle(slmField);
    sptField = sum(input_intensity.*exp(1i*(-slmPh+sptPhOffset)),[1 2])/length(xj(:));
    wt = mean(abs(sptField))/abs(sptField).*wt;
    
    sptPhOffset = angle(sptField);%./sptField;
    slmField = sum(exp(1i*(-sptPhOffset-sptPh)).*wt,3);
    
    
    
end

phase = gather(slmPh);
%%
% phase = reshape(phase(slm_mask),1152/2,1920/2);

% phase = imresize(phase,[1152 1920],'nearest');

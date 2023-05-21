load images
figure windowstyle docked
imshow(image160626191120,[]), colorbar
figure windowstyle docked
imshow(image160626192006,[]), colorbar
figure windowstyle docked
imshow(lockin160626192006,[]), colorbar
sb = image160626191120;

%%
ghsb = fdag(cat(3,-lockin160626192006,zeros(300,300)),mean(sb(:))); % fourier domain antigradient
asb = antigradient(cat(3,-lockin160626192006,zeros(300,300)),mean(sb(:))); % farneback antigradient

%%
figure windowstyle docked, imshow(ghsb,[]), colorbar
figure windowstyle docked, imshow(asb,[]), colorbar

%%
% ysb = [zeros(1,size(sb,2)); diff(image160626191120,1,1)];
% xsb = [zeros(size(sb,1),1), diff(image160626191120,1,2)];
[xsb ysb] = gradient(-image160626191120);
figure windowstyle docked, imshow(ysb,[]), colorbar

%%
ghsb = fdag(cat(3,-ysb,zeros(300,300)),mean(sb(:))); % fourier domain antigradient
asb = antigradient(cat(3,-ysb,zeros(300,300)),mean(sb(:))); % farneback antigradient

%%
figure windowstyle docked, imshow(ghsb,[]), colorbar
figure windowstyle docked, imshow(asb,[]), colorbar

% figure windowstyle docked, imshow([sb-ghsb],[]), colorbar
% sqrt(sum((sb(:)-ghsb(:)).^2))

% figure windowstyle docked, imshow([sb-asb],[]), colorbar
% sqrt(sum((sb(:)-asb(:)).^2))

% figure windowstyle docked, imshow([sb-ghsb sb-asb],[]), colorbar
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb = image160626191120;

fsb = ifftshift(fft2(fftshift(sb)));
fsb = real(fsb);
figure windowstyle docked
imshow(real(fsb),[]), colorbar
subplot 211, plot(sb(end/2,:))
subplot 212, plot(fsb(end/2,:))
ysb = [zeros(1,size(sb,2)); diff(sb,1,1)];
xsb = [zeros(size(sb,1),1), diff(sb,1,2)];
[xsb ysb] = gradient(sb);
clf, imshow(xsb,[])
plot(xsb(end/2,:))
ghsb = fdag(cat(3,ysb,xsb),mean(sb(:))); % fourier domain antigradient
asb = antigradient(cat(3,ysb,xsb),mean(sb(:))); % farneback antigradient

%%
figure windowstyle docked, imshow(sb,[]), colorbar
figure windowstyle docked, imshow(ghsb,[]), colorbar
figure windowstyle docked, imshow(asb,[]), colorbar

figure windowstyle docked, imshow([sb-ghsb],[]), colorbar
sqrt(sum((sb(:)-ghsb(:)).^2))

figure windowstyle docked, imshow([sb-asb],[]), colorbar
sqrt(sum((sb(:)-asb(:)).^2))

figure windowstyle docked, imshow([sb-ghsb sb-asb],[]), colorbar


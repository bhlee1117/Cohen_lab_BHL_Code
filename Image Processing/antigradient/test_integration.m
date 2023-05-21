sb = poissrnd(abs(phantom)*1000)+100;
sb = sb + randn(size(sb))*10;
sb = abs(sb);
figure, imshow(sb,[])

whos sb
fsb = ifftshift(fft2(fftshift(sb)));
fsb = real(fsb);
imshow(real(fsb),[]), colorbar
subplot 211, plot(sb(end/2,:))
subplot 212, plot(fsb(end/2,:))
xsb = [zeros(1,size(sb,2)); diff(sb,1,1)];
ysb = [zeros(size(sb,1),1), diff(sb,1,2)];
[xsb ysb] = gradient(sb);
imshow(xsb,[])
plot(xsb(end/2,:))
ghsb = fdag(cat(3,ysb,xsb),mean(sb(:))); % fourier domain antigradient
asb = antigradient(cat(3,ysb,xsb),mean(sb(:))); % farneback antigradient

%%
figure, imshow(sb,[]), colorbar
figure, imshow(ghsb,[]), colorbar
figure, imshow(asb,[]), colorbar

figure, imshow([sb-ghsb],[]), colorbar
sqrt(sum((sb(:)-ghsb(:)).^2))

figure, imshow([sb-asb],[]), colorbar
sqrt(sum((sb(:)-asb(:)).^2))

figure, imshow([sb-ghsb sb-asb],[]), colorbar


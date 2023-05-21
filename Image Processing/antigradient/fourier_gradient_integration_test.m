synthbump = load('sampledbumps','h2p5ir1p');
fncell = fieldnames(synthbump);
synthbump = eval(['synthbump.' fncell{1} ';']);
pad = 200;
subs = 2;
sb = synthbump([2*ones(1,pad) 1:end-1 end*ones(1,pad)],[2*ones(1,pad) 1:end-1 end*ones(1,pad)]);
sb = sb(ceil(1:subs:(end*3/4)),1:subs:end)/100;
% sb = randn(480,640);
% imshow(sb,[])
% whos sb
% fsb = ifftshift(fft2(fftshift(sb)));
% fsb = real(fsb);
% imshow(real(fsb),[]), colorbar
% subplot 211, plot(sb(end/2,:))
% subplot 212, plot(fsb(end/2,:))
% xsb = [zeros(1,size(sb,2)); diff(sb,1,1)];
% ysb = [zeros(size(sb,1),1), diff(sb,1,2)];
[xsb ysb] = gradient(sb);
% imshow(xsb,[])4
% plot(xsb(end/2,:))
ghsb = fdag(cat(3,ysb,xsb),mean(sb(:)));
asb = antigradient(cat(3,ysb,xsb),mean(sb(:)));

%%
figure, imshow(sb,[]), colorbar
figure, imshow(ghsb,[]), colorbar
figure, imshow(asb,[]), colorbar

figure, imshow([sb-ghsb],[]), colorbar
sqrt(sum((sb(:)-ghsb(:)).^2))

figure, imshow([sb-asb],[]), colorbar
sqrt(sum((sb(:)-asb(:)).^2))

figure, imshow([sb-ghsb sb-asb],[]), colorbar


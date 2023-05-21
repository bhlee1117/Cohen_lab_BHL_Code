h=get(gca,'children');p = h.CData;
% figure(20);clf;imagesc(p);colorbar
p_n=mat2gray(p);figure(20);imagesc(p_n);colorbar;axis equal
function mov_pca_recon=pca_ica_bhl(mov_bin_s)
tic
n_pc = 16;
%mov_bin_s=mov_bin_s(15:end-15,15:end-15,1:3000);
%mov_xcorr=toimg(tovec(mov_bin_s(:,:,:)).*tovec(mean(mov,3)),124,152);
%covMat = tovec(mov_bin_s(:,:,1:3000))*tovec(mov_bin_s(:,:,1:3000))';
 %[V_pc, D_pc,U_pc] = eigs(double(covMat),n_pc);
%in = gpuArray(double(tovec(mov_bin_s)));
[V_pc, D_pc,U_pc] = svds(double(tovec(mov_bin_s)),n_pc);
%[V_pc, D_pc,U_pc] = svds(in,n_pc);
clear in
V_pc = gather(V_pc);
D_pc = gather(D_pc);
U_pc = gather(U_pc);
 toc
D_pc = diag(D_pc); 
% D = D(end:-1:1);
% V = V(:,end:-1:1);
vSign = sign(max(V_pc) - max(-V_pc));   % make the largest value always positive
V_pc = V_pc.*vSign;
toc

idx_pc = 1:16;
% pc_sub = (V_pc(:,idx_pc)'*tovec(mov_bin_s))'; %PCA
% pc_sub = pc_sub-mean(pc_sub); %PCA
pc_sub = U_pc; %svds
dt = 1.27e-3;
% figure(9);clf
% plot((1:2900)*dt,pc_sub./range(pc_sub(501:end-500,:))+(1:length(idx_pc)))
% ylim([0 length(idx_pc)+1])

figure(10);clf
for ii = 1:16
    subplot(4,4,ii)
    imshow2(toimg(V_pc(:,idx_pc(ii)),size(mov_bin_s,1),size(mov_bin_s,2)),[prctile(V_pc(:,idx_pc(ii)),1) prctile(V_pc(:,idx_pc(ii)),99)])
    title(num2str(ii))
end


idx_pc = input('Select PC');
[icasig, A,W] = fastica((V_pc(:,idx_pc)'*tovec(mov_bin_s)));
n_ic = size(W,1);
V_ic = (W * V_pc(:,idx_pc)')';
ic_sub = (V_ic(:,1:n_ic)'*tovec(mov_bin_s))';
mov_pca_recon = SeeResiduals(mov_bin_s,ic_sub);
end
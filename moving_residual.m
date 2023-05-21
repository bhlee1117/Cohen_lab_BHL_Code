function res_trace=moving_residual(traces,mcTrace,pca_trace,window)

seg=[1:window/2:size(traces,2)-window]'; wind=[0:window];
ns=size(traces,1);
for j=1:size(seg,1)
for jj=1:size(traces,1)
c=corrcoef(traces(jj,seg(j)+wind),mcTrace(seg(j)+wind,1));
Corrs(jj,j,1)=c(2,1);
c=corrcoef(traces(jj,seg(j)+wind), mcTrace(seg(j)+wind,2));
Corrs(jj,j,2)=c(2,1);
c=corrcoef(traces(jj,seg(j)+wind),pca_trace(seg(j)+wind));
Corrs(jj,j,3)=c(2,1);
end
end
Corrs=abs(Corrs); 
res_trace=traces;
pca_hi=pca_trace-movmean(pca_trace,50);

for j=find(sum(Corrs(:,:,1)>0.8,1)>0 | sum(Corrs(:,:,2)>0.8,1)>0 | sum(Corrs(:,:,3)>0.8,1)>0)
    res_trace(:,seg(j)+wind)=squeeze(SeeResiduals(reshape(traces(:,seg(j)+wind),ns,1,[]),mcTrace(seg(j)+wind,:)'));
    res_trace(:,seg(j)+wind)=squeeze(SeeResiduals(reshape(res_trace(:,seg(j)+wind),ns,1,[]),mcTrace(seg(j)+wind,:).^2'));
    res_trace(:,seg(j)+wind)=squeeze(SeeResiduals(reshape(res_trace(:,seg(j)+wind),ns,1,[]),(mcTrace(seg(j)+wind,1).*mcTrace(seg(j)+wind,2))'));
    %res_trace(:,seg(j)+wind)=squeeze(SeeResiduals(reshape(res_trace(:,seg(j)+wind),ns,1,[]),pca_trace(seg(j)+wind)'));
    %res_trace(:,seg(j)+wind)=squeeze(SeeResiduals(reshape(res_trace(:,seg(j)+wind),ns,1,[]),pca_hi(seg(j)+wind)'));
    
  %res_trace(:,seg(j,1):seg(j,2))=squeeze(SeeResiduals(reshape(Result{i}.traces_res(:,seg(j,1):seg(j,2)),135,1,[]),pca_trace(seg(j,1):seg(j,2))'));
end
end
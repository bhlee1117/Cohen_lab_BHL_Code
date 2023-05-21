function roi=dmd_pat2roi(dmd_mask_sequence_rois,sz)
roi=zeros(sz(1),sz(2));
t=find(isnan(dmd_mask_sequence_rois(:,1)));
for i=1:length(t)
   if t(1)==6
roi=roi+roipoly(roi,dmd_mask_sequence_rois(t(i)-5:t(i)-1,1),dmd_mask_sequence_rois(t(i)-5:t(i)-1,2));
   else
       try
roi=roi+roipoly(roi,dmd_mask_sequence_rois(t(i-1)+1:t(i)-1,1),dmd_mask_sequence_rois(t(i-1)+1:t(i)-1,2));                  
       catch
roi=roi+roipoly(roi,dmd_mask_sequence_rois(1:t(i)-1,1),dmd_mask_sequence_rois(1:t(i)-1,2));       
       end
   end
end
roi=roi>0;
end
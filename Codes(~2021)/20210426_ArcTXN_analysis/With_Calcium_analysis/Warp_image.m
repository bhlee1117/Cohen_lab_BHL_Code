function Warp_im=Warp_image(ImageToWarp,tform2)
for i=1:size(ImageToWarp,3)
[transVis xdata ydata] = imtransform(ImageToWarp(:,:,i), tform2);
transVisWithOffset = imtransform(ImageToWarp(:,:,i), tform2, 'XData', [1 (size(ImageToWarp(:,:,i),2)+tform2.tdata.T(3,1))],'YData', [1 (size(ImageToWarp(:,:,i),1)+ tform2.tdata.T(3,2))]);
transVisBuffrd = uint16(zeros(size(ImageToWarp(:,:,i))));
transVisBuffrd(1:size(transVisWithOffset,1),1:size(transVisWithOffset,2)) = transVisWithOffset;
Transformation_Matrix=tform2.tdata.T(1:2,1:2);
Translation_Vector=tform2.tdata.T(3,1:2)';
transVisBuffrd = transVisBuffrd((1:size(ImageToWarp,1)), (1:size(ImageToWarp,2)));
Warp_im(:,:,i)=transVisBuffrd;
end
end

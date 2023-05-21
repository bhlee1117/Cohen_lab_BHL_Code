function out = img_process(filename1, filename2, firstimage, lastimage)
% function out = img_process(filename, firstimage, lastimage)
% Displays the difference between the image arrays
% Plots the time trajectory of the intensities
% Returns the g factor; g=(I1-I2)/[(I1+I2)/2]
% Malcolm Campbell 7/14/08
img1=read_tifs2(filename1,firstimage,lastimage);
img2=read_tifs2(filename2,firstimage,lastimage);
avgimg1=mean(img1,3);
avgimg2=mean(img2,3);
diff=avgimg1-avgimg2;
intens1=squeeze(mean(mean(img1,1),2));
intens2=squeeze(mean(mean(img2,1),2));
t=firstimage+1:lastimage+1;
imtool(diff)
plot(t,intens1,t,intens2)
out=mean((intens1-intens2)./((intens1+intens2)/2));
m1 = vm('R:\S1\FOV1\s1');
m2 = vm('R:\S1\FOV1\s2');
m3 = vm('R:\S1\FOV1\s3');
m4 = vm('R:\S1\FOV1\s4');
m5 = vm('R:\S1\FOV1\s5');
m = ((m1+m2+m3)/3+m4+m5)/3;
% s1 - s3 same FOV
% s1 - s5 using zero zernike correction
% s6 using -20 defocus, -40 spherical
% s7 using widefield patterned illumination
m6 = vm('R:\S1\FOV1\s6');
m7 = vm('R:\S1\FOV1\s7');

md = permute(m7(:,:,end:-1:1).data,[3 1 2]);
md = imresize(md,[6*1.6/2.25 1].*[m.frames m.cols]);
figure
moviefixsc(md)
figure
moviefixsc(permute(md,[1 3 2]))

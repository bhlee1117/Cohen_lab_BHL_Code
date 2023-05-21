function volt = AOTF_LUT(perct)

lut = load('AOTF_LUT.txt','sheet1','A:A');
[lut_perct,lut_idx] = unique(lut(:,2));
lut_v = lut(lut_idx,1);
volt = interp1(lut_perct,lut_v,perct);
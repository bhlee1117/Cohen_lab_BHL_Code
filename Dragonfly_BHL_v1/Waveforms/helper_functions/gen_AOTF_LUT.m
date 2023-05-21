function gen_AOTF_LUT(sheet,col1,col2)
[f,p] = uigetfile('*.xlsx');
cur_dir = pwd;
cd(p)
lut_v = xlsread(f,sheet,[col1 ':' col1]);
lut_perct = xlsread(f,sheet,[col2 ':' col2]);
lut = [lut_v lut_perct];
save('AOTF_LUT.txt','lut','-ascii');
cd(cur_dir)
pat = load('D:\Code\Dragonfly_yq_v1\slm_control\1024x1024.mat','pic_fft_out');
pat = pat.pic_fft_out;

hSLM.project(pat)
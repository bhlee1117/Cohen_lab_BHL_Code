function S=deconvolve_BH(trace,rising_time,decay_time,noise)

tau_d = decay_time/30;  
tau_r = rising_time/30; 
nMax = 100; 
pars = [tau_d, tau_r]; 
kernel = create_kernel('exp2', pars, nMax); 

[chat, S, kernel_fit, iters] = deconvCa(trace, kernel, 2, true, false,noise);

end
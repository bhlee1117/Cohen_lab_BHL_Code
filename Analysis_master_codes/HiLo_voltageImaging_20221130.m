%HiLo_voltage imaging

clear
[fname fpath] = uigetfile('*.*');
mov=single(vm(fpath));
ref_im=movmax(mov,3,3);
for i=1:50%size(mov,3)
H=hilospeckle(ref_im(:,:,i),mov(:,:,i));
H.targetThickness=3;        
H.runHiLo;
mov_HiLo(:,:,i)=H.HiloFinal;
end
   
%%
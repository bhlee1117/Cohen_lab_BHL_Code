%% Center to wide pattern blue intensity calibration

V = 4; % LED voltage used for calibration test

fp = 5; % camera exposure frame per pattern


[fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
mov = vm(fpath);
[~,intens] = clicky_faster(mov,mov(1:5).mean);

f_c = mean(intens(2:fp))-100;
f_w = mean(intens(fp+2:fp*2))-100;

V_w = f_c/f_w * V
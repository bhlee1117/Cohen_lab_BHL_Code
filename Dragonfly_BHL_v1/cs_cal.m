function V_center=cs_cal(V_wide)
[fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
mov = vm(fpath);
[~,F]=clicky_autobg(mov,mov(1:end/2).mean);


% FA1/FA2=1/QSa*Ic/Is
Ic = power_cal(V_wide,'blue','control');
Is = Ic;
Qsa=Ic/Is*F(end/2+1)/F(end/2-1);
V_center = power_cal(Qsa*power_cal(V_wide,'blue','control'),'blue','intensity');
end


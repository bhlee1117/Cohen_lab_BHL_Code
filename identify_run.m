function R=identify_run(Arduino_data,thres)
dif=Arduino_data(2:end,2)-Arduino_data(1:end-1,2);
dif(dif>9000)=0;
T=cumsum(dif);
Velocity=T(2:end)-T(1:end-1);
R=Velocity>thres;
end
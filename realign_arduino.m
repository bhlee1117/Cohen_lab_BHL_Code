function Arduino_align=realign_arduino(Arduinos)
for i=1:length(Arduinos)
Ard_data=Arduinos{i};
reward_spot=mean(Ard_data(Ard_data(:,3)==1,2));
lap_length(i)=max(Ard_data(:,2))-min(Ard_data(:,2));
Ard_data(:,2)=Ard_data(:,2)-reward_spot;
Ard_data(find(Ard_data(:,2)<-lap_length(i)/2),2)=lap_length(i)+Ard_data(find(Ard_data(:,2)<-lap_length(i)/2),2);
Ard_data(find(Ard_data(:,2)>lap_length(i)/2),2)=Ard_data(find(Ard_data(:,2)>lap_length(i)/2),2)-lap_length(i);
Arduino_align{i}=Ard_data;
end
end





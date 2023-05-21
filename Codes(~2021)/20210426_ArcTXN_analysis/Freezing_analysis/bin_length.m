function output=bin_length(input_data)

a=input_data(1,1);
r=1;
s=1;
output=zeros(1,2);
for k=2:size(input_data,1) %length of track

if input_data(k,1)==a
    s=s+1;
else
    output(r,a+1)=s;
    r=r+1;
    s=1;
    a=input_data(k,1);
end
end 

output(r,a+1)=s;
    
end
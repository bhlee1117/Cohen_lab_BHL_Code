function new_labels=switchlabel(original_labels)

unique_labels = unique(original_labels);
% Initialize a new label vector
new_labels = zeros(size(original_labels));

% Loop through the unique values and assign increasing integers
for i = unique_labels'
    newind_t(i)=find(original_labels == unique_labels(i),1);
end    
newind_t=sort(newind_t,'ascend');

for i = 1:length(unique_labels)
  new_labels(original_labels == original_labels(newind_t(i))) = i;
end
end
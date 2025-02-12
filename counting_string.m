function ordinal_string=counting_string(num_list,prefix)
if nargin<2
    prefix=[];
end
  % Initialize an empty string
    ordinal_string = [];
    
    % Iterate through the list of numbers
    for i = 1:length(num_list)
        num = num_list(i);
        
        % Determine the correct ordinal suffix
        if num >= 11 && num <= 13
            suffix = 'th';
        else
            switch mod(num, 10)
                case 1
                    suffix = 'st';
                case 2
                    suffix = 'nd';
                case 3
                    suffix = 'rd';
                otherwise
                    suffix = 'th';
            end
        end
        
        % Add the number with its suffix to the ordinal_string
        ordinal_string{i} = strcat(sprintf('%s%d%s', prefix, num, suffix));
        
    end
end
function device = load_device_data(dd, device_identifier, id_type)
    if ~exist('id_type', 'var')
        id_type = "Device_Name";
    end
    for i = 1:length(dd)
       if isfield(dd{i}, id_type) && eval(sprintf("dd{i}.%s", id_type)) == device_identifier
           device = dd{i};
           return;
       end
    end
    device = {};
end
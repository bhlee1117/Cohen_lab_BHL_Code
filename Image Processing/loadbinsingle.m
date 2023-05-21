        function obj = loadbinsingle(fname,sz)
            if numel(sz) == 1
                sz = [sz sz];
            end
            fid = fopen(fname,'r');
            data = fread(fid,'*float32');
            data = reshape(data,sz(1),sz(2),[]);
            fclose(fid);
            obj = vm(data);
        end

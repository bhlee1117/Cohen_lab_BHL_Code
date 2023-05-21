%ffzpad     FFT friendly pad or crop.
% 
%   2016 Vicente Parot
%   Cohen Lab - Harvard University
%
        function m = ffzpad(m,ts)
            % resize canvas (pad or crop) to centered target size ts = [rows cols]
            m = permute(ffzpadrows(m,ts(1)),[2 1 3:ndims(m)]);
            m = permute(ffzpadrows(m,ts(2)),[2 1 3:ndims(m)]);
            function m = ffzpadrows(m,rows)
                [sr, sc, sf] = size(m);
                dr = rows - sr;
                if dr < 0
                    m = m(1-floor(dr/2):end+ceil(dr/2),:,:);
                elseif dr > 0
                    m = [zeros(ceil(dr/2),sc,sf); m; zeros(floor(dr/2),sc,sf)];
                end    
            end
        end

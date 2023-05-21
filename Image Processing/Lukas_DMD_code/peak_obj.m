classdef peak_obj < handle
    properties
        n %number of pulses
        dti %length of each pulse
        v % one over frequency
        length %length of peak in #of frames
        loc %location of the peak
        data
    end
    
    methods 
        function peak_obj = peak_obj(location, length, n, dti)
            peak_obj.loc = location; 
            peak_obj.length = length; 
            peak_obj.v = floor(peak_obj.length / n); 
            peak_obj.n = n; 
            peak_obj.dti = dti; 
            
            if (peak_obj.dti * peak_obj.n > peak_obj.length)
                error('peak:inArgs', 'Problem: mismatch between number of peaks, lengths of peaks, and length of waveform'); 
            end
            
            peak_obj.data = zeros(1, peak_obj.length); 
            for i = 1:n
                offset = peak_obj.v; 
                indices_begin = (offset*(i-1) + 1);
                indices_end = indices_begin + (peak_obj.dti - 1); 
                peak_obj.data(1, indices_begin:indices_end) = ones(1, peak_obj.dti);
                peak_obj.data = boolean(peak_obj.data); 
            end
        end
    end
end
        
        
        
            
        
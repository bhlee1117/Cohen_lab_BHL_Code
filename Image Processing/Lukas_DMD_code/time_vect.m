classdef time_vect < handle
    properties
        n %number of peaks
        peaks %vector containing peaks objects
        length %length in frames
        data
    end
    methods 
        function time_vect = time_vect(length, varargin)
            time_vect.n = nargin - 1;
            time_vect.length = length; 
            
            if nargin > 0
                for i = 1:time_vect.n; 
                    time_vect.add_peak(varargin{i}); 
                end
            end
            
            time_vect.data = boolean(zeros(1,length));
        end
        
        function add_peak(time_vect, varargin)
            
            for i = 1:size(varargin)
                time_vect.peaks{i} = varargin{i}; 
            end
            time_vect.n = time_vect.n + size(varargin); 
        end
        
        function value = get.data(obj)
            for i = 1:size(obj.peaks)
                obj_indices = obj.peaks{i}.loc:(obj.peaks{i}.loc + obj.peaks{i}.length - 1); 
                if max(obj_indices > obj.length)
                    error('time_vect:peakLengths', 'You have peaks that are outside bounds of time_vect object')
                end
                obj.data(1, obj_indices) = obj.peaks{i}.data;
            end
            value = obj.data; 
        end
    end
end
            
        
        
              
                
    
        
        
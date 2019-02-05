classdef Block
    % class for the description of one block, with position of different
    % plattens
    
    % the origin of the block coordinates is at the West-South-Bottom corner (WSB)
    
    properties
    
         L_E = 250
         L_N = 250
         L_T = 250
         sizes = [250, 250, 250]; % in mm
         
         faces_offset;        
    end
    
    methods
    
        % constructor
        function obj = Block(sizes)
        
            % size - vector of length 
            obj.sizes = sizes;
           
        end
        
        % methods -> plot block ?
        
        
    end
    
    
end

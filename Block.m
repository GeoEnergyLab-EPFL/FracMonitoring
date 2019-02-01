classdef Block
    % class for the description of one block, with position of different
    % plattens
    
    % the origin of the block coordinates is at the West-South-Bottom corner (WSB)
    
    properties
    
         L_E = 250
         L_N = 250
         L_T = 250
         sizes = [250, 250, 250]; % in mm
         
         faces(2,6) cell % mapping between faces and platten object -> like 2d array of characters
         
         faces_offset;
          
    end
    
    methods
    
        % constructor
        function obj = Block(faces,sizes)
        
            % size - vector of length 
            obj.sizes = sizes;
            obj.faces = faces;
           
        end
        
        % methods -> plot block ?
        
        
    end
    
    
end

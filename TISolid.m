classdef TISolid < Solid
    % description
    %
    % classs for Transversely isotropic solid properties
    %
    % todo : Put the computation of group velocity as function of ray angle
    %        as methods of this class, also add a C_ij get function
    %  
    properties
        
        c11;
        c33;
        c13;
        c55;
        c66;
    end
    
    properties (Dependent)
        
        c44;
        
        C_ij;%
        
    end
    
    methods
        
        % constructor
        function obj=TISolid(density,c11,c33,c13,c55,c66,varargin)
            
            obj=obj@Solid(density,varargin);
            
            obj.c11=c11;
            obj.c33=c33;
            obj.c13=c13;
            obj.c55=c55;
            obj.c66=c66;
            
        end
        
        
        % Group velocity should move here !
        
        
        % code up C_ij matrix get
        
    end
end

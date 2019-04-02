classdef Solid
    %  description:
    %  Mother class for solid properties
    %  contains density and solid name (optional)
    %
    % Brice Lecampion - 2019
    %
    %
    properties
        
        name="blank";
        density=1000; % kg/m^3
        
    end
    
    methods
        
        % constructor
        function obj=Solid(density,varargin)
            
            obj.density=density;
            
            if (isempty(varargin{1}))
                obj.name="blank";
            else
                obj.name=varargin{1}{1};
            end
            
        end
        
    end
end


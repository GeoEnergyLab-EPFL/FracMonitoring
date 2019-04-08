classdef Fluid
    % description:
    %
    % a class for fluids used in HF experiments, containing the physical
    % properties of the fluid relevant for the injection
    
    properties
        
        rho;
        compressibility;
        Vp;
        viscosity;
        
    end
    
    methods
        % constructor
        function obj = Fluid(fluidtype)
            if strcmp(fluidtype,'glycerol')
                % rho and Vp from JASA paper (doi:10.1121/1.1907292)
                obj.rho = 1260;% density glycerol
                obj.Vp = 1960;  % velocity glycerol
                obj.compressibility = 0;
                obj.viscosity = 0; % in Pa.s
            elseif strcmp(fluidtype,'siliconeA')
                obj.rho = 1050;%
                obj.Vp = 1350;  %
                obj.compressibility = 0;
                obj.viscosity = 0; % in Pa.s.
            else
                fprintf('\nError: unknown fluid!\n');
                return
            end
        end
        % end constructor
    end
end
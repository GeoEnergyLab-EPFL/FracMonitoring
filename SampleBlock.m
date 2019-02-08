classdef SampleBlock
    % class for the description of one block, with position of different
    % plattens
    
    % the origin of the block coordinates is at the West-South-Bottom corner (WSB)
    
    properties
         L_E double % length in E-W direction in mm
         L_N double % length in N-S direction in mm
         L_T double % length in T-B direction in mm
         sizes = [250, 250, 250]; % vector with the three lengths above
         faces_offset % offsets with respect to WSB origin
    end
    
    methods
        % constructor
        function obj = SampleBlock(sizes)
            % size - vector of length 
            obj.sizes = sizes;
            obj.L_E=sizes(1);
            obj.L_N=sizes(2);
            obj.L_T=sizes(3);
        end
        
        % methods 
        
        % plot block geometry in 3D
        function fig_handle = blockplot3D(obj,varargin)
            % open figure from passed handle if it exists
            if ~isempty(varargin)
                if isgraphics(varargin{1})
                    fig_handle = figure(varargin{1});
                else
                    fig_handle = figure;
                end
            else
                fig_handle = figure;
            end
            hold on
            % plot block edges
            plot_cuboid(obj.sizes,[0 0 0],fig_handle)
            axis equal tight
        end
        
    end
    
end

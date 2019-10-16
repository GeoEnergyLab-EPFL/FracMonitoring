classdef SampleBlock
    % description
    % class encapsulating data of one block smaple, with position of different
    % plattens
    %
    % the origin of the block coordinates is at the West-South-Bottom corner (WSB)
    
    properties
        L_E double % length in E-W direction in m
        L_N double % length in N-S direction in m
        L_T double % length in T-B direction in m
        sizes = [0.25, 0.25, 0.25]; % vector with the three lengths above
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
        
        % METHODS
        % plot block geometry in 3D
        % todo : add label N, S , E , W , T, B  for each sides
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
        
        % PLOT BLOCK WITH PLATTENS
        function fig_handle = plotblockwithplattens(obj,platten_list,varargin)
            
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
            
            for i=1:length(platten_list)
                fig_handle=plattenplot3D(platten_list(i),obj,fig_handle);
            end
            
        end
        
        % Dong Liu -- 16/09/2019
        % PLOT ONE FACE OF THE BLOCK WITH CORRESPONDING TRANSDUCERS
        function fig_handle = plotside2D(obj,sidemarker,varargin)
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
            switch sidemarker
                case 'S'
                    l_x = obj.L_N;
                    l_y = obj.L_T;
                    xlabel('Easting');
                    ylabel('Height');
                case 'N'
                    l_x = obj.L_N;
                    l_y = obj.L_T;
                    xlabel('Easting'); 
                    ylabel('Height');
                case 'E'
                    l_x = obj.L_E;
                    l_y = obj.L_T;
                    xlabel('Northing');
                    ylabel('Height');
                case 'W'
                    l_x = obj.L_E;
                    l_y = obj.L_T;
                    xlabel('Northing');
                    ylabel('Height');
                case 'T'
                    l_x = obj.L_N;
                    l_y = obj.L_E;
                    xlabel('Easting');
                    ylabel('Northing');
                case 'B'
                    l_x = obj.L_N;
                    l_y = obj.L_E;
                    xlabel('Easting');
                    ylabel('Northing');
                otherwise
                    disp('Wrong input for direction indicator, N,S,E,W,T,B are accepted as input');
                    return;
            end
            hold on
            title([sidemarker ' platten'],'FontSize',30);
            hold on
            rectangle('Position',[0 0 l_x l_y]);
            axis equal tight
            
        end
        
        % PLOT ONE FACE OF THE BLOCK WITH CORRESPONDING TRANSDUCERS
        % Dong Liu -- 16/09/2019
        function fig_handle = plotside2Dwithplattens(obj,platten_list,sidemarker,varargin)
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
            fig_handle = plotside2D(obj,sidemarker,fig_handle);
            hold on
            for i=1:length(platten_list)
                if (platten_list(i).face == sidemarker)
                    fig_handle = plattenplot2Dside(platten_list(i),sidemarker,fig_handle);
                end
            end

        end
        
        
    end
    
end

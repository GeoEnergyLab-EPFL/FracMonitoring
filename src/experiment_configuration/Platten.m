classdef Platten
    % description:
    % A class for the definition of platten and their position during an experiment
    % This class defines the platten object, that includes info on the
    % physical platten and hole positions as well as how the platten is
    % used in the global context of the reaction frame
    
    properties
        type char       % cross or grid
        id(1,1) char    % id of the platten
        xy_holes(:,2) double     % positions of the transducer holes in local coordinates
        face(1,1) char
        
        ep_x double % vector coordinates defining local e_x in the global block system of coordinates
        ep_y double % vector coordinates defining local e_y in the global block system of coordinates
        ep_z double % outward normal in the global block system of coordinates (pointing out of the block)
        offset
    end
    
    properties (Dependent)
        n_holes % number of holes in platten
        R % rotation matrix from local to global coord system
    end
    
    methods
        % constructor
        function obj = Platten(ID,face,block_data) % ,varargin
            
            obj.id = ID; % platten id
            if ismember(face,{'N','S','E','W','T','B'})
                obj.face = face;
            else
                disp('Wrong sample face for platten')
            end
            
            switch ID
                case {'A','B'}
                    % it's a cross
                    obj.type='cross';
                    
                    % positions of holes on platen A and B in mm
                    dcross = 18;    % top platen cross spacing
                    ddiag1 = 34;    % top platen first diag spacing
                    ddiag2 = 28;    % top platen additional diag spacing
                    % compute positions in patten frame of reference,
                    % origin is center of platten
                    x_tmp = [zeros(1,10), [5 4 3 2 1]*dcross,...
                        -[1 2 3 4 5]*dcross, [-(ddiag1+2*ddiag2)...
                        -(ddiag1+ddiag2) -ddiag1 ddiag1 ddiag1+ddiag2,...
                        ddiag1+2*ddiag2 -(ddiag1+2*ddiag2) -(ddiag1+ddiag2)...
                        -ddiag1 ddiag1, ddiag1+ddiag2 ddiag1+2*ddiag2]/sqrt(2)]';
                    y_tmp = [[5 4 3 2 1]*dcross, -[1 2 3 4 5]*dcross,...
                        zeros(1,10), [(ddiag1+2*ddiag2) (ddiag1+ddiag2)...
                        ddiag1 ddiag1 ddiag1+ddiag2, ddiag1+2*ddiag2...
                        -(ddiag1+2*ddiag2) -(ddiag1+ddiag2) -ddiag1 -ddiag1,...
                        -(ddiag1+ddiag2) -(ddiag1+2*ddiag2)]/sqrt(2)]';
                    % write to property
                    obj.xy_holes = [x_tmp, y_tmp];
                    
                case {'C','D','E','F'}
                    % it's a grid
                    obj.type='grid';
                    
                    % positions of holes on platen C, D, E and F in mm
                    dhoriz = 50;    % side platen horiz spacing
                    dvert1 = 20;    % side platen first vert spacing
                    dvert2 = 30;    % side platen farther vert spacing
                    % compute positions
                    x_tmp = repmat([-dhoriz 0 dhoriz],1,4)';
                    y_tmp = [ones(1,3)*(dvert1+dvert2), ones(1,3)*dvert1,...
                        ones(1,3)*-dvert1, ones(1,3)*-(dvert1+dvert2)]';
                    % write to property
                    obj.xy_holes = [x_tmp, y_tmp];
            end
            
            %             % add varargin
            %             nVarargs = length(varargin);
            %             if (nVarargs==2)
            %                 obj.offset_x = varargin{1};
            %                 obj.offset_y = varargin{2};
            %             end
            
            switch (obj.face)
                case 'N'
                    obj.ep_z = [0,1,0];
                    obj.ep_x = [-1,0,0];
                    obj.ep_y = [0,0,1]; % local coord vect in global coord system
                    obj.offset = [block_data.sizes(1)/2., block_data.sizes(2),...
                        block_data.sizes(3)/2.];
                    
                case 'S'
                    
                    obj.ep_z = [0,-1,0];
                    obj.ep_x=[1,0,0];
                    obj.ep_y=[0,0,1];
                    obj.offset = [block_data.sizes(1)/2.,0., block_data.sizes(3)/2.];
                    
                case 'E'
                    
                    obj.ep_z = [1,0,0];
                    obj.ep_x=[0,1,0];
                    obj.ep_y=[0,0,1];
                    
                    obj.offset = [block_data.sizes(1) ,block_data.sizes(2)/2.,...
                        block_data.sizes(3)/2.];
                    
                case 'W'
                    obj.ep_z = [-1,0,0];
                    obj.ep_x = [0,-1,0];
                    obj.ep_y = [0,0,1];
                    obj.offset = [0. ,block_data.sizes(2)/2., block_data.sizes(3)/2.];
                    
                case 'T'
                    obj.ep_z = [0,0,1];
                    obj.ep_x = [1,0,0];
                    obj.ep_y = [0,1,0];
                    obj.offset = [block_data.sizes(1)/2. ,block_data.sizes(2)/2.,...
                        block_data.sizes(3)];
                    
                case 'B'
                    obj.ep_z = [0,0,-1];
                    obj.ep_x = [-1,0,0];
                    obj.ep_y = [0,1,0];
                    obj.offset = [block_data.sizes(1)/2. ,block_data.sizes(2)/2.,...
                        0.];
            end
            
        end
        % end constructor
        
        
        % METHODS
        % method for dependant property "n_holes" number of holes
        function value = get.n_holes(obj)
            value = length(obj.xy_holes);
        end
        
        % method for dependant property "R" rotation matrix
        function value = get.R(obj)
            % create rotation
            value = [obj.ep_x',obj.ep_y',obj.ep_z'];
        end
        
        % plot platten geometry in 2D
        function fig_handle = plattenplot2D(obj,varargin)
            % open figure from passed handle if it exists
            if ~isempty(varargin)
                if isgraphics(varargin{1})
                    fig_handle = figure(varargin{1});
                    hold on
                else
                    fig_handle = figure;
                end
            else
                fig_handle = figure;
            end
            plot(obj.xy_holes(:,1),obj.xy_holes(:,2),'ok')
            hold on
            % add hole numbering
            holetxt = cellstr(num2str((0:(obj.n_holes-1))'));
            txtoffset = 2;
            text(obj.xy_holes(:,1)+txtoffset,obj.xy_holes(:,2),holetxt)
            
            % add edges of platten
            plot([-125 -125],[-125 125],'k')
            plot([-125 125],[125 125],'k')
            plot([125 125],[-125 125],'k')
            plot([-125 125],[-125 -125],'k')
            axis equal tight
        end
        
        % plot platten geometry in 3D 
        function fig_handle = plattenplot3D(obj,block_data,varargin)
            % compute 3D locations
            xyz_loc = (obj.R*[obj.xy_holes(:,1), obj.xy_holes(:,2),...
                zeros(size(obj.xy_holes(:,1)))]')'+obj.offset;
            
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
            % plot location of transducer holes for platten object
            plot3(xyz_loc(:,1),xyz_loc(:,2),xyz_loc(:,3),'ok')
            % add hole numbering
            holetxt = cellstr(num2str((0:(obj.n_holes-1))'));
            %txtoffset = 2;
            text(xyz_loc(:,1),xyz_loc(:,2),xyz_loc(:,3),holetxt)
            % add block edges ---- this will be dublicated plenty of times
            % ;(
            plot_cuboid(block_data.sizes,[0 0 0],fig_handle)
            axis equal tight
            % add platten letter ID in top-left corner of outer platten
            % face (as is physically stamped)
            letterIDpos = (obj.R*[-100 100 0]')'+obj.offset;
            text(letterIDpos(1),letterIDpos(2),letterIDpos(3),obj.id)
        end
        
    end
    
end


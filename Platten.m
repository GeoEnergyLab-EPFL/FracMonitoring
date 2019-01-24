classdef Platten
    % A class for the definition of platten and their position during an experiment
    %   Detailed explanation goes here
    
    properties
        type ; % cross or grid
        id ;   % id of the block
        
        xy_holes;   % positions of the transducer holes in local coordinates
    end
    
    
    properties (Dependent)
        n_holes;
    end
    
    methods
        
        % constructor
        function obj = Platten(ID) % ,varargin
            
            obj.id = ID; % platten id
            
            switch ID
                case {'A','B'}
                    % it's a cross
                    obj.type='cross';
                    
                    % positions of holes on platen A and B in mm
                    dcross = 18;    % top platen cross spacing
                    ddiag1 = 34;    % top platen first diag spacing
                    ddiag2 = 28;    % top platen additional diag spacing
                    % compute positions
                    x_tmp = [zeros(1,10), [5 4 3 2 1]*dcross, -[1 2 3 4 5]*dcross,...
                        [-(ddiag1+2*ddiag2) -(ddiag1+ddiag2) -ddiag1 ddiag1 ddiag1+ddiag2,...
                        ddiag1+2*ddiag2 -(ddiag1+2*ddiag2) -(ddiag1+ddiag2) -ddiag1 ddiag1,...
                        ddiag1+ddiag2 ddiag1+2*ddiag2]/sqrt(2)]';
                    y_tmp = [[5 4 3 2 1]*dcross, -[1 2 3 4 5]*dcross, zeros(1,10),...
                        [(ddiag1+2*ddiag2) (ddiag1+ddiag2) ddiag1 ddiag1 ddiag1+ddiag2,...
                        ddiag1+2*ddiag2 -(ddiag1+2*ddiag2) -(ddiag1+ddiag2) -ddiag1 -ddiag1,...
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
                    x_tmp = repmat([-dhoriz 0 dhoriz],1,4);
                    y_tmp = [ones(1,3)*(dvert1+dvert2), ones(1,3)*dvert1,...
                        ones(1,3)*-dvert1, ones(1,3)*-(dvert1+dvert2)];
                    % write to property
                    obj.xy_holes = [x_tmp, y_tmp];
            end
            
            % get number of objects from length of coordinate property
            obj.n_holes = length(obj.xy_holes);
            
            %             % add varargin
            %             nVarargs = length(varargin);
            %             if (nVarargs==2)
            %                 obj.offset_x=varargin{1};
            %                 obj.offset_y=varargin{2};
            %             end
            
            
        end
        % end constructor
        
        function value=get.n_holes(obj)
            value = size(obj.xy_holes,1);
        end
        
        % METHODS
        % 1 method -> plot platten geometry
        
    end
    
end


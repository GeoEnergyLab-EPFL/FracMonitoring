classdef Transducers
    
    properties
        serial(:,1) int32   % vector containing transducer serial numbers
        type(:,1) char      % vector containing either 'S' or 'R'
        channel(:,1) int32  % vector containing the transducer channel
        platten(:,1) char   % vector containing the platten ID on which the transducer is
        local_id(:,1) int32 % vector containing the hole ID of the transducer
        orientation(:,1) double   % vector containing the orientation of the transducer
        % cable for Shear transducer, angle in radians with respect to the 
        % platten local coordinate system
    end
    
    properties (Dependent)
        n_transducers
        n_sources
        n_receivers
    end
    
    
    methods
        
        % constructor
        function obj = Transducers(serial,type,channel,platten,loc_id,orient)
            % put checks on all vectors length !!!
            if ~isequal(length(serial),length(type),length(channel),...
                    length(platten),length(loc_id),length(orient))
                fprintf('\nError: properties should all have the same length!\n');
                return
            else
                obj.serial = serial;
                obj.type = type;
                obj.channel = channel;
                obj.platten = platten;
                obj.local_id = loc_id;
                obj.orientation = orient;
            end
            
        end
        
        % methods for dependant properties
        % number of transducers
        function val = get.n_transducers(obj)
            val = length(obj.serial);
        end
        
        % number of sources
        function val = get.n_sources(obj)
            channelchar = char(obj.channel);
            val = sum(ismember(channelchar(:,1),'S'));
        end
        
        % number of receivers
        function val = get.n_receivers(obj)
            channelchar = char(obj.channel);
            val = sum(ismember(channelchar(:,1),'R'));
        end
        
        % METHODS
        % calculate global coordinates (in block reference) of each transducer
        function xyzTransd = calc_global_coord(obj,platten_list)
            % platten_list: list of platten objects - max length 6
            % block_data: a Block object  containing the size of the block
            % sizes: L_E, L_N, L_T  which means that the origin of the
            % block coordinates is at the corner (WSB)
            %
            % loop on number of piezo ...
            xyzTransd = zeros(obj.n_transducers,3);
            
            for ii = 1:obj.n_transducers
                % find corresponding platten
                p = 1;
                while p <= length(platten_list)
                    if strcmp(obj.platten(ii),platten_list{p}.id)
                        % get xy in local platten coordinates system
                        xyloc = platten_list{p}.xy_holes(obj.local_id(ii)+1,:); % vector of length 2
                        
                        % ONLY AFTER OFFSET ADDED IN PLATTEN CLASS
                        %                         % add platten offset
                        %                         xyloc(1) = xyloc(1)+platten_list{p}.offset_x;
                        %                         xyloc(2) = xyloc(2)+platten_list{p}.offset_y;
                        
                        % local coordinates
                        x_l = xyloc(1);
                        y_l = xyloc(2);
                        z_l = 0.;
                       
                        % in global coordinates
                        xyzTransd(ii,:)= ((platten_list{p}.R)*[x_l;y_l;z_l])' + platten_list{p}.offset ;
                        break
                    end
                    p=p+1;
                    
                end
                
                % error if platten not found in list
                if (p>length(platten_list))
                    disp('transducers location NOT found !! - inconsistent data');
                    xyzTransd(ii,:)= [-999.,-999.,-999.];
                end
                
            end
            
        end
        
        % 3D plot of the transducer locations
        function fig_handle = transducerplot3D(obj,platten_list,varargin)
            % open figure from passed handle if it exists
            % optional argument 1 is figure handle
            % optional argument 2 is plotting style
            if ~isempty(varargin)
                narg = length(varargin);
                if isgraphics(varargin{1})
                    fig_handle = figure(varargin{1});
                else
                    fig_handle = figure;
                end
            else
                fig_handle = figure;
            end
            hold on
            xyzTransd = calc_global_coord(obj,platten_list);
            plotstyleS = 'r.';
            plotstypeR = 'b.'
            if narg>=2
                if ischar(varargin{2})
                    plotstyleS = varargin{2};
                    plotstyleR = varargin{2};
                end
            end
            plot3(xyzTransd(:,1),xyzTransd(:,2),xyzTransd(:,3),plotstyle)
        end
        
        % distances for all source-receiver pairs
        function dists = transducerdists(obj,platten_list)
            xyzTransd = calc_global_coord(obj,platten_list);
            dx = xyzTransd(:,1)-xyzTransd(:,1)';
            dy = xyzTransd(:,2)-xyzTransd(:,2)';
            dz = xyzTransd(:,3)-xyzTransd(:,3)';
            dists = sqrt(dx.^2+dy.^2+dz.^2);
        end
        
        % directions of rays for all source-receiver pairs
        function directions = transducerdirections(obj,platten_list)
            xyzTransd = calc_global_coord(obj,platten_list);
            dx = xyzTransd(:,1)-xyzTransd(:,1)'; % diff in x direction
            dy = xyzTransd(:,2)-xyzTransd(:,2)'; % diff in y direction
            dz = xyzTransd(:,3)-xyzTransd(:,3)'; % diff in z direction
            dirtmp = cat(3,dx,dy,dz); % array with 3D differences
            directions = dirtmp./vecnorm(A,2,3); % normalized
        end
        
    end
    
end


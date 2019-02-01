classdef Transducers
    % Transudcers class goes here
    %   
    
    properties
        
        serial(:,1) int32   % vector containing transducer serial numbers
        type(:,1) char      % vector containing either 'P' or 'S'
        channel(:,1) string % vector string, how the transducers are used (for example "S24")
        platten(:,1) char   % vector containing the platten ID on which the transducer is   
        local_id(:,1) int32 % vector containing the hole ID of the transducer
        % add orientation for S transducers???

    end
    
    properties (Dependent)
    
        n_transducers
        n_sources
        n_receivers

    end
    
    
    methods
        
        % constructor
        function obj = Transducers(serial,type,channel,platten,loc_id)
            % put checks on all vectors length !!!
            if ~isequal(length(serial),length(type),length(channel),...
                    length(platten),length(loc_id))
                fprintf('\nError: properties should all have the same length!\n');
                return
            else
            obj.serial = serial;
            obj.type = type;
            obj.channel = channel;
            obj.platten = platten;
            obj.local_id = loc_id;
            end

        end
        
        % methods for dependant properties
        function val = get.n_transducers(obj)
            val = length(obj.serial);  
        end
        
        function val = get.n_sources(obj)
            channelchar = char(obj.channel);
            val = sum(ismember(channelchar(:,1),'S'));
        end
        
        function val = get.n_receivers(obj)
            channelchar = char(obj.channel);
            val = sum(ismember(channelchar(:,1),'R'));
        end
        
        % calculate global coordinates (in block reference) of each transducer
        function xyz = calc_global_coord(obj,platten_list,block_data)
            % platten_list: list of platten objects - max length 6
            % block_data: a Block object  containing the size of the block
            % sizes: L_E, L_N, L_T  which means that the origin of the
            % block coordinates is at the corner (WSB)
            % 
            % loop on number of piezo ...
            xyz = zeros(obj.n_transducers,3);
            
            for i = 1:obj.n_transducers
                
                % find corresponding platten 
                p = 1;
                while p <= length(platten_list)
                    if (obj.platten(i) == platten_list{p}.id)
                        
                        % get xy in local platten coordinates system 
                        xyloc = platten_list{p}.xy_holes(obj.local_id(i),:); % vector of length 2
                        disp(xyloc);
                        
                        % CHANGE TO FACE FROM BLOCK OBJECT
                        disp(platten_list{p}.face);
                        
% ONLY AFTER OFFSET ADDED IN PLATTEN CLASS                        
%                         % add platten offset
%                         xyloc(1) = xyloc(1)+platten_list{p}.offset_x;
%                         xyloc(2) = xyloc(2)+platten_list{p}.offset_y;
                        
                        % local coordinates
                        x_l = xyloc(1);
                        y_l = xyloc(2);
                        z_l = 0.;
                        % create global 
                        % 1 get face normal
                        
                            % add offset from platten position.... 
                        % be careful platten local system is w.r. platten
                        % center
                        disp([' platten face: ', platten_list{p}.face]);
                        switch (platten_list{p}.face)  
                            case 'N'
                                n=[0,1,0];
                                offset = [block_data.sizes(1)/2., block_data.sizes(2), block_data.sizes(3)/2.];
                            case 'S'
                                n=[0,-1,0];
                                 offset = [block_data.sizes(1)/2.,0., block_data.sizes(3)/2.];
                            case 'E'
                                n=[1,0,0];
                                 offset = [block_data.sizes(1) ,block_data.sizes(2)/2., block_data.sizes(3)/2.];
                            case 'W'
                                n=[-1,0,0];
                                 offset = [0. ,block_data.sizes(2)/2., block_data.sizes(3)/2.];
                            case 'T'
                                n=[0,0,1];
                                offset = [block_data.sizes(1)/2. ,block_data.sizes(2)/2., block_data.sizes(3)];
                            case 'B'
                                n=[0,0,-1];
                                  offset = [block_data.sizes(1)/2. ,block_data.sizes(2)/2., 0.];
                        end 
                        
                        % create rotation
                        R = [platten_list{p}.xloc',  platten_list{p}.yloc', n'];
                       
                        % in global aux variables x_g, y_g and z_g
                        xyz(i,:)= (R*[x_l;y_l;z_l])' + offset ;
                        
                    end
                    p=p+1;
                    
                end
                
                % error if platten not found in list
                if (p>length(platten_list))
                    disp('transducers location NOT found !! - inconsistent data');
                    xyz(i,:)= [-999.,-999.,-999.];
                end
   
            end
      
        end
        
%         % transducer plotting method
%         function fighandle = plot_transducer(obj,platten_list,block_data)
%             xyz = calc_global_coord(obj,platten_list,block_data);
%             
%         end

    end    
    
end


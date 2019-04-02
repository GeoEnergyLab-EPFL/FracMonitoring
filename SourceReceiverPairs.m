classdef SourceReceiverPairs
    % description: 
    % class encapsulating data associated with source - receiver pairs 
    %
    % Brice Lecampion - 2019
    %
    %
    
    properties
        
        SRmap;  % matrix of integer containing the mapping   Source Receiver  size( n_pair by 2)
        
        XS_XR ;  % matrix of real containing coord of source , cord of receivers - n_pair rows by 6 colum
        
        wave_type;% vector of length n_pair with type of wave - string ?
        
    end
    
    properties (Dependent)
        
        distances;
        n_pairs;
        
    end
    
    methods
        
        % constructor
        function obj=SourceReceiverPairs(TransducerObj,platten_list,my_map,varargin)
            %note we put wave_type in varargin - as sometimes one may not
            %need it - in that case we just put dummy values
            
            [~,nc] = size(my_map);
            
            if (nc~=2)
                disp('error in mapping input - not 2 columns!!');
                return
            end
            
            % check on sources in the map
            C = intersect(TransducerObj.channel(1:TransducerObj.n_sources),unique(my_map(:,1)));
            if ( length(C)<length(unique(my_map(:,1))) )
                disp('error - some source in the map are not registered in the transducers array' );
                return
            end
            % check on receivers in the map
            
            C = intersect(TransducerObj.channel(1:TransducerObj.n_sources),unique(my_map(:,2)));
            if ( length(C)<length(unique(my_map(:,2))) )
                disp('error - some receivers in the map are not registered in the transducers array' );
                return
            end
            
            obj.SRmap= my_map;
            
            % 1 compute the location of all the transducers - even if we
            % do not need them all
            xyzTransd = calc_global_coord(TransducerObj,platten_list);
            
            xyz_source = xyzTransd(my_map(:,1),:);
            xyz_receiver = xyzTransd(my_map(:,2)+TransducerObj.n_sources,:);
            
            obj.XS_XR=[xyz_source xyz_receiver];
            
            if (length(varargin)==1)
                [nr,nc]=size(varargin{1});
                if (nc==1) && (nr==length(my_map(:,1)) )
                    obj.wave_type=varargin{1};
                else
                    disp('error in wave type numbered entered - not consistent with number of S-R pairs in the map entered');
                    return;
                end
                
            end
            
        end
        
        
        % method to get distances between all R-S pairs
        function distances=get.distances(obj)
            
            df=obj.XS_XR(:,1:3)-obj.XS_XR(:,4:6);
            distances=sqrt(df(:,1).^2+df(:,2).^2+df(:,3).^2);
            
        end
        
        % number of source - receiver pairs in the object
        function n_pairs = get.n_pairs(obj)
             
            [n_pairs ~ ] =size(obj.SRmap);
            
            
        end
        
        % get distance for a given pair in the map
        function d=getDistancePairI(obj,i)
            
            df=obj.XS_XR(i,1:3)-obj.XS_XR(i,4:6);
            d=norm(df);
            
        end
        
        % function to get distance for a range of pairs in the map
        function d=getDistancePairRange(obj,r)
            
            df=obj.XS_XR(r,1:3)-obj.XS_XR(r,4:6);
            d=sqrt(df(:,1).^2+df(:,2).^2+df(:,3).^2);
            
        end
        
        
        % function to plot the rays of all pairs in the map
        function fig_handle=plotdirectrays(obj,varargin)
          
           % narg = length(varargin);
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
           
            
            % Here do a arrow ?
            % WTF seems obvious but hard to do in matlab
            
            % by default Source  POSITION in RED - Receiver in Blue 
            % and the ray in black
 
            for i=1:obj.n_pairs
                plot3(obj.XS_XR(i,[1 4])',obj.XS_XR(i,[2 5])',obj.XS_XR(i,[3 6])','k.-.','linewidth',2)
                hold on;
            end
            plot3(obj.XS_XR(:,1),obj.XS_XR(:,2),obj.XS_XR(:,3),'r.','MarkerSize',18);
            hold on
            plot3(obj.XS_XR(:,4),obj.XS_XR(:,5),obj.XS_XR(:,6),'b.','MarkerSize',18);
            
%             
            axis equal
        
        end
        
    end
    
end



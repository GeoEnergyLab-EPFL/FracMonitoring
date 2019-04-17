classdef Transducers
    % description
    %
    % This class encapsulates all data of an array of transducers:
    % - serial #
    % - transducer type: Source or Receiver
    % - channel id (integer) -> giving corresponding index of the data
    % matrix of a sequence
    % - platten id : platten where the transducer is located
    % - local_id : hole  # where the transducer is located on the platten
    % - orientation: orientation of calble (for shear wave transducer)
    % - wave_mode: type of transducer either longitudinal or shear
    %
    %
    properties
        serial(:,1) int32   % vector containing transducer serial numbers
        type(:,1) char      % vector containing either 'S' or 'R'
        channel(:,1) int32  % vector containing the transducer channel
        platten(:,1) char   % vector containing the platten ID on which the transducer is
        local_id(:,1) int32 % vector containing the hole ID of the transducer
        orientation(:,1) double   % vector containing the orientation of the transducer
        % cable for Shear transducer, angle in radians with respect to the
        % platten local coordinate system
        wave_mode ; % vector containing the type of transducers either Longitudinal (0) or Shear (1)
        
    end
    
    properties (Dependent)
        
        n_transducers;
        n_sources;
        n_receivers;
        
    end
    
    methods
        
        % constructor
        function obj = Transducers(serial,type,channel,platten,loc_id,orient)
            % put checks on all vectors length !!!
            if ~isequal(length(serial),length(type),length(channel),length(platten),length(loc_id),length(orient))
                fprintf('\nError: properties should all have the same length!\n');
                return
            else
                
                %%% REORDERING IN SEQUENCE Source  0-31 / Receiver 0-31
                % instead of by Platten
                
                [~, iaux] = sort(type); % order by alphabetical so first R then S
                typechar = char(type);
                nr = sum(ismember(typechar(:,1),'R'));
                ns = length(type)-nr;
                iaux_s = iaux(nr+1:nr+ns);
                iaux_r = iaux(1:nr);
                
                [~, is] = sort(channel(iaux_s));
                [~, ir] = sort(channel(iaux_r));
                
                re_order = [iaux_s(is); iaux_r(ir)];
                
                % check on unicity of channel for the sources
                if (length(unique(channel(iaux_s)))~=ns)
                    disp('error some duplicate channel in sources ');
                    return
                end
                % check on unicity of channel for the receivers
                if (length(unique(channel(iaux_r)))~=nr)
                    disp('error some duplicate channel in receivers ');
                    return
                end
                
                obj.channel = channel(re_order);
                obj.serial = serial(re_order);
                obj.type = type(re_order);
                obj.platten = platten(re_order);
                obj.orientation = orient(re_order);
                obj.local_id = loc_id(re_order);
                
                obj.wave_mode =(obj.orientation~=0);
                % to be checked with  the serial numbers ?
                % a priori it seems to be ok - prototype P is 1608142 and
                % proto shear is 1608143 after all shear transducers are above 1609XXXX
                
            end
            
        end
        
        % methods for dependant properties
        % number of transducers
        function val = get.n_transducers(obj)
            val = length(obj.serial);
        end
        
        % number of sources
        function val = get.n_sources(obj)
            val = length(find(obj.type=='S'));
        end
        
        % number of receivers
        function val = get.n_receivers(obj)
            val = length(find(obj.type=='R'));
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
                    
                    if strcmp(obj.platten(ii),platten_list(p).id)
                        % get xy in local platten coordinates system
                        xyloc = platten_list(p).xy_holes(obj.local_id(ii)+1,:); % vector of length 2
                        
                        % ONLY AFTER OFFSET ADDED IN PLATTEN CLASS
                        %                         % add platten offset
                        %                         xyloc(1) = xyloc(1)+platten_list{p}.offset_x;
                        %                         xyloc(2) = xyloc(2)+platten_list{p}.offset_y;
                        
                        % local coordinates
                        x_l = xyloc(1);
                        y_l = xyloc(2);
                        z_l = 0.;
                        
                        % in global coordinates
                        xyzTransd(ii,:)= ((platten_list(p).R)*[x_l;y_l;z_l])' + platten_list(p).offset ;
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
        
        % method returning the  transducers channel of transducers
        % located in a given platten either source or receiver
        % this function is sort of depreciated now with the SourceReceiver
        % Pairs object
        function channels_on_platten = Transducers_on_platten(obj,platten_char,s_r_char)
            
            % check that platten_char is either  on obj
            pl_obj=unique(obj.platten);
            
            if isempty(intersect(platten_char,pl_obj))
                disp('The letter of the platten given is not in the transducer objet ! error ');
                return
            end
            
            kc=find(obj.platten==platten_char);
            
            if (s_r_char ~='S') && (s_r_char ~='R')
                disp('The letter given for source or receiver is neither S or R ! error ');
                return
            end
            
            if (s_r_char=='S')
                channels_on_platten=kc(obj.type(kc)==s_r_char);
            else
                channels_on_platten=kc(obj.type(kc)==s_r_char)-obj.n_sources; % receiver case offset id by number of sources
            end
            
        end
        
        % 3D plot of the transducer locations
        function fig_handle = transducerplot3D(obj,platten_list,varargin)
            % open figure from passed handle if it exists
            % optional argument 1 is figure handle
            % optional argument 2 is plotting style
            
            narg = length(varargin);
            
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
            
            xyzTransd = calc_global_coord(obj,platten_list);
            plotstyleS = 'r.';
            plotstyleR = 'b.';
            if narg>=2
                if ischar(varargin{2})
                    plotstyleS = varargin{2};
                    plotstyleR = varargin{2};
                end
                
            end
            
            plot3(xyzTransd(1:obj.n_sources,1),xyzTransd(1:obj.n_sources,2),...
                xyzTransd(1:obj.n_sources,3),plotstyleS);
            plot3(xyzTransd(obj.n_sources+1:end,1),xyzTransd(obj.n_sources+1:end,2),...
                xyzTransd(obj.n_sources+1:end,3),plotstyleR);
        end
        
        % distances for all source-receiver pairs
                % this function should be a method of source -receiver pair
                % object -> unneeded here
        function dists = transducerdists(obj,platten_list)
            xyzTransd = calc_global_coord(obj,platten_list);
            dx = xyzTransd(:,1)-xyzTransd(:,1)';
            dy = xyzTransd(:,2)-xyzTransd(:,2)';
            dz = xyzTransd(:,3)-xyzTransd(:,3)';
            dists = sqrt(dx.^2+dy.^2+dz.^2);
        end
        
        % directions of rays for all source-receiver pairs
        % this function should be a method of source -receiver pair object
        function directions = transducerdirections(obj,platten_list)
            xyzTransd = calc_global_coord(obj,platten_list);
            dx = xyzTransd(:,1)-xyzTransd(:,1)'; % diff in x direction
            dy = xyzTransd(:,2)-xyzTransd(:,2)'; % diff in y direction
            dz = xyzTransd(:,3)-xyzTransd(:,3)'; % diff in z direction
            dirtmp = cat(3,dx,dy,dz); % array with 3D differences
            directions = dirtmp./vecnorm(A,2,3); % normalized
        end
        
        
        % Construct  S-R pairs object without the intra-platten pairs
        function objpair=AllexceptintraPairs(TransducerObj,platten_list)
            % - do not choose the S-R which are on the same platten  only the S-R pairs
            % All S-R pairs beside the one on the same platten than the source
            n_s = length(find(TransducerObj.type=='S'));
            n_r = length(find(TransducerObj.type=='R'));
            myMap=zeros((n_s*n_r),2);
            wave_type_mat=zeros((n_s*n_r),2);
            ch=linspace(1,n_r,n_r)';
            k=1;
            for i=1:n_s
                % myTransducers are ordered Source - Then receiver
                source_platten=TransducerObj.platten(i);
                R_on_platten=Transducers_on_platten(TransducerObj,source_platten,'R');
                
                other_r=setdiff(ch,R_on_platten);
                
                myMap(k:k+length(other_r)-1,1)=i;
                myMap(k:k+length(other_r)-1,2)=other_r;
                wave_type_mat(k:k+length(other_r)-1,1)= TransducerObj.wave_mode(i);
                wave_type_mat(k:k+length(other_r)-1,2)=TransducerObj.wave_mode(other_r);
                
                k=k+length(other_r);
                
            end
            
            
            myMap=myMap(1:k-1,:);
             wave_type=[];
            
            for i=1:k-1   
                if (~wave_type_mat(i,1))
                    wt='P';             
                else
                    wt='S';
                end
                
                if (~wave_type_mat(i,2))
                    wave_type=[wave_type ; char(strcat(wt,'P')) ];
                else
                    wave_type=[wave_type ; char(strcat(wt,'S')) ];
                end
            end
            
            objpair=SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        function objpair=TwoPlattenPairs(TransducerObj,platten_list,p_1,p_2)
            
            % create SRPairs object for opposite platten taking the source
            % of platten p_1 and the receivers on platten p_2
            
            % checks on string p_1 and p_2
            n_p=length(platten_list);
            
            pl_f=char(zeros(n_p,1));
            for p=1:n_p
                pl_f(p)= platten_list(p).face;
            end
            
            if isempty(intersect(pl_f,p_1)) || isempty(intersect(pl_f,p_2))
                disp(' error in given platten face or corresponding platten list input!');
                return
            end
            
            p_first = find(pl_f==p_1) ;
            p_two = find(pl_f==p_2) ;
            
            % find source on platten
            ks_1= find(((TransducerObj.platten==platten_list(p_first).id) .* (TransducerObj.type=='S'))==1);
            % kr_1= find(((TransducerObj.platten==platten_list(p_first).id) .* (TransducerObj.type=='R'))==1)-TransducerObj.n_sources;
            
            s_type=TransducerObj.wave_mode(ks_1);
            % find receiver on opposite platten
            %  ks_2= find(((TransducerObj.platten==platten_list(p_two).id) .* (TransducerObj.type=='S'))==1);
            kr_2= find(((TransducerObj.platten==platten_list(p_two).id) .* (TransducerObj.type=='R'))==1) -TransducerObj.n_sources;
            r_type=TransducerObj.wave_mode(kr_2);
            myMap=zeros((length(ks_1) )*(length(kr_2)),2);
            %            myMap=zeros((length(ks_1)+length(ks_2))*(length(kr_1)+length(kr_2)),2);
            wave_type_mat=zeros((length(ks_1) )*(length(kr_2)),2);
            
            k=1;
            for i=1:length(ks_1)
                myMap(k:k+length(kr_2)-1,1)=ks_1(i);
                myMap(k:k+length(kr_2)-1,2)=kr_2;
                wave_type_mat(k:k+length(kr_2)-1,1)=s_type(i);
                wave_type_mat(k:k+length(kr_2)-1,2)=r_type;
                k=k+length(kr_2);
            end
            
            %             for i=1:length(ks_2)
            %                 myMap(k:k+length(kr_1)-1,1)=ks_2(i);
            %                 myMap(k:k+length(kr_1)-1,2)=kr_1;
            %                 k=k+length(kr_1);
            %             end
            myMap=myMap(1:k-1,:); 
            wave_type=[];
            
            for i=1:k-1   
                if (~wave_type_mat(i,1))
                    wt='P';             
                else
                    wt='S';
                end
                
                if (~wave_type_mat(i,2))
                    wave_type=[wave_type ; char(strcat(wt,'P')) ];
                else
                    wave_type=[wave_type ; char(strcat(wt,'S')) ];
                end
            end
            
            objpair=SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        
        % constructor of  S-R pairs for opposite platten pairs
        function objpair = AllOppositePairs(TransducerObj,platten_list)
            % - do not choose the S-R which are on the same platten  only the S-R pairs
            % All S-R pairs beside the one on the same platten than the source
            n_s = length(find(TransducerObj.type=='S'));
            n_r = length(find(TransducerObj.type=='R'));
            myMap = zeros((n_s*n_r),2);
            wave_type_mat=zeros((n_s*n_r),2);
            n_p = length(platten_list);
            pl_f = char(zeros(n_p,1));
            for p = 1:n_p
                pl_f(p) = platten_list(p).face;
            end
            
            k = 1;
            for p = 1:n_p
                
                % find opposite platten
                
                switch platten_list(p).face
                    case 'B' % find 'T'
                        p_opp = find(pl_f=='T') ;
                        
                    case 'T' % find 'B'
                        p_opp = find(pl_f=='B') ;
                    case 'E'
                        p_opp = find(pl_f=='W') ;
                    case 'W'
                        p_opp = find(pl_f=='E') ;
                    case 'N'
                        p_opp = find(pl_f=='S') ;
                    case 'S'
                        p_opp = find(pl_f=='N') ;
                    otherwise
                        disp('error in platten face and opposite- check your inputs')
                end
                
                % find source on platten
                ks= find(((TransducerObj.platten==platten_list(p).id) .* (TransducerObj.type=='S'))==1);
                s_type=TransducerObj.wave_mode(ks);
                % find receiver on opposite platten
                kr= find(((TransducerObj.platten==platten_list(p_opp).id) .* (TransducerObj.type=='R'))==1) -TransducerObj.n_sources;
                 r_type=TransducerObj.wave_mode(kr);
                % loop on sources on that platten & add in map
                for i=1:length(ks)
                    myMap(k:k+length(kr)-1,1)=ks(i);
                    myMap(k:k+length(kr)-1,2)=kr;
                     wave_type_mat(k:k+length(kr)-1,1)=s_type(i);
                    wave_type_mat(k:k+length(kr)-1,2)=r_type;
                    k=k+length(kr);
                end

            end
            
            myMap=myMap(1:k-1,:);
            
            wave_type=[];
            
            for i=1:k-1   
                if (~wave_type_mat(i,1))
                    wt='P';             
                else
                    wt='S';
                end
                
                if (~wave_type_mat(i,2))
                    wave_type=[wave_type ; char(strcat(wt,'P')) ];
                else
                    wave_type=[wave_type ; char(strcat(wt,'S')) ];
                end
            end
            
            
            objpair=SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
    end
    
end

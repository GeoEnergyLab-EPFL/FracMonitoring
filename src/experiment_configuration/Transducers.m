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
            % optional argument 3 is channel numbering (matlab style from
            % 1), boolean 
            
            narg = length(varargin);
            
            if ~isempty(varargin)&&~isempty(varargin{1})&&isgraphics(varargin{1})
                fig_handle = figure(varargin{1});
            else
                fig_handle = figure;
                axis square
                xlabel('Easting (m)')
                ylabel('Northing (m)')
                zlabel('Height (m)')
            end
            hold on
            
            xyzTransd = calc_global_coord(obj,platten_list);
            plotstyleS = 'ro';
            plotstyleR = 'bo';
            % selected plot style if option selected
            if narg>=2
                if ~isempty(varargin{1})&&ischar(varargin{2})
                    plotstyleS = varargin{2};
                    plotstyleR = varargin{2};
                end
                
            end
            % change marker size
            mkrsize = 10;
            % plot sources in 3D
            plot3(xyzTransd(1:obj.n_sources,1),xyzTransd(1:obj.n_sources,2),...
                xyzTransd(1:obj.n_sources,3),plotstyleS,'MarkerSize',mkrsize);
            % plot receivers in 3D
            plot3(xyzTransd(obj.n_sources+1:end,1),xyzTransd(obj.n_sources+1:end,2),...
                xyzTransd(obj.n_sources+1:end,3),plotstyleR,'MarkerSize',mkrsize);
            % add black marker to identify shear transducers
            plot3(xyzTransd(obj.wave_mode==1,1),xyzTransd(obj.wave_mode==1,2),...
                xyzTransd(obj.wave_mode==1,3),'.k')
            
            % add transducer numbering if option selected
            if narg>=3
                if ~isempty(varargin{3})&&isnumeric(varargin{3})
                    offset = 0.004;     % offset
                    text(xyzTransd(:,1)+offset,xyzTransd(:,2)+offset,xyzTransd(:,3)+offset,...
                        num2cell(obj.channel));
                end
            end
            
        end
        
        % 2D plot of the transducer locations for one face of the block
        % Dong Liu -- 16/09/2019
        function fig_handle = transducerplot2D(obj,platten_list,sidemarker,varargin)
            % open figure from passed handle if it exists
            % optional argument 1 is figure handle
            % optional argument 2 is plotting style
            % optional argument 3 is channel numbering (matlab style from
            % 1), boolean 
            
            narg = length(varargin);
            
            if ~isempty(varargin)&&~isempty(varargin{1})&&isgraphics(varargin{1})
                fig_handle = figure(varargin{1});
            else
                fig_handle = figure;
            end
            hold on
            
            xyzTransd = calc_global_coord(obj,platten_list);
            plotstyleS = 'ro';
            plotstyleR = 'bo';
            % selected plot style if option selected
            if narg>=2
                if ~isempty(varargin{1})&&ischar(varargin{2})
                    plotstyleS = varargin{2};
                    plotstyleR = varargin{2};
                end
                
            end
            % change marker size
            mkrsize = 10;
            
            switch sidemarker
                case 'N'
                    i=1;
                    j=3;
                case 'S'
                    i=1;
                    j=3;
                case 'E'
                    i=2;
                    j=3;
                case 'W'
                    i=2;
                    j=3;
                case 'T'
                    i=1;
                    j=2;
                case 'B'
                    i=1;
                    j=2;
                otherwise
                    disp('Wrong input for side indicator')
                    return;
            end
            
            % add transducer numbering if option selected
            if narg>=3
                if ~isempty(varargin{3})&&isnumeric(varargin{3})
                    offset = 0.004; % offset
                    for i_platten =1:length(platten_list)
                        if platten_list(i_platten).face == sidemarker
                            i_plattenid=platten_list(i_platten).id;
                        end
                    end
                    for i_transducer = 1:obj.n_transducers
                        if obj.platten(i_transducer) == i_plattenid
                            % plot sources and receivers in 2D
                            if i_transducer<=obj.n_sources
                                plotstyle=plotstyleS;
                            else
                                plotstyle=plotstyleR;
                            end
                            plot(xyzTransd(i_transducer,i),...
                            xyzTransd(i_transducer,j),plotstyle,'MarkerSize',mkrsize);   
                            
                            % add black marker to identify shear transducers
                            if obj.wave_mode(i_transducer)==1
                                plot(xyzTransd(i_transducer,i),...
                                xyzTransd(i_transducer,j),'.k')
                            end
                            
                            text(xyzTransd(i_transducer,i)+offset,...
                                xyzTransd(i_transducer,j)+offset,...
                        num2cell(obj.channel(i_transducer)));
                        end

                    end
                end
            end
            
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
        function objpair = AllexceptintraPairs(TransducerObj,platten_list)
            % - do not choose the S-R which are on the same platten  only the S-R pairs
            % All S-R pairs beside the one on the same platten than the source
            n_s = length(find(TransducerObj.type=='S'));
            n_r = length(find(TransducerObj.type=='R'));
            myMap = zeros((n_s*n_r),2);
            wave_type_mat = zeros((n_s*n_r),2);
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
            
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        function objpair = TwoPlattenPairs(TransducerObj,platten_list,p_1,p_2)
            
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
            
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        % constructor of one source or receiver on one platten and all the 
        % receivers or sources on the opposite platten
        % Dong Liu -- 17/09/2019
        function objpair = SinglePlattenPairs(TransducerObj,platten_list,p_1,p_2,t_1,SRindicator)
            
            % create SRPairs object for opposite platten taking the source
            % of platten p_1 and the receivers on platten p_2
            
            % checks on string p_1 and p_2
            TwoplattenPairs=TwoPlattenPairs(TransducerObj,platten_list,p_1,p_2);
            if SRindicator == 'S'
                idx=find(TwoplattenPairs.SRmap(:,1) == t_1);
            elseif SRindicator == 'R'
                idx=find(TwoplattenPairs.SRmap(:,2) == t_1);
            else
                disp('Wrong input, only S or R is acceptable');
            end
            myMap = TwoplattenPairs.SRmap(idx,:);
            wave_type = TwoplattenPairs.wave_type(idx,:);
            
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        
        % constructor of  S-R pairs for opposite platten pairs
        function objpair = AllPairsOppositePlattens(TransducerObj,platten_list)
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
            
            
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        % constructor of S-R pairs from two plattens
        % S from platten_1 with R from platten_2 and S from platten_2 with
        % R from platten_1
        % Dong Liu -- 04/11/2019
        function objpair = TwoPlattenPairsAll(TransducerObj,platten_list,p_1,p_2)
            
            pair_1 = TwoPlattenPairs(TransducerObj,platten_list,p_1,p_2);
            pair_2 = TwoPlattenPairs(TransducerObj,platten_list,p_2,p_1);
            myMap = [pair_1.SRmap; pair_2.SRmap];
            wave_type = [pair_1.wave_type;pair_2.wave_type];
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        % constructor of S-R pairs from the same platten
        % Dong Liu -- 04/11/2019
        function objpair = OnePlattenPairsAll(TransducerObj,platten_list,p_1)
            
            % create SRPairs object for one platten
            objpair = TwoPlattenPairs(TransducerObj,platten_list,p_1,p_1);    
        end        
        
        % constructor of S-R pairs for the diffraction with one side chosen
        % and all the other transducers from the top and the bottom
        % plattens
        % Dong Liu -- 04/11/2019
        function objpair = SidePlattenPairsAll(TransducerObj,platten_list,p_1)
            
            % create SRPairs object for neighboring platten taking the source
            % of platten p_1 and the receivers on platten top and bottom
            % and taking the receivers of platten p_1 and sources on the
            % platten top and bottom
            pair_1 = TwoPlattenPairsAll(TransducerObj,platten_list,p_1,'T');
            pair_2 = TwoPlattenPairsAll(TransducerObj,platten_list,p_1,'B');
            myMap = [pair_1.SRmap; pair_2.SRmap];
            wave_type = [pair_1.wave_type;pair_2.wave_type];
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);   
        end
        

        % constructor of  S-R pairs for the same platten pairs in order to
        % look at the reflection data for transducers on the top and the
        % bottom plattens
        % DongLiu -- 31/10/2019
        % we can also look at the pairs for the two faced plattens in other
        % two directions Dong Liu -- 01/11/2019
        % TODO: add the direction option
        function objpair = AllPairsSamePlattens(TransducerObj,platten_list)
            % choose the S-R which are on the same platten
            
            n_sB = length(find((TransducerObj.type=='S').*(TransducerObj.platten=='A')));
            n_rB = length(find((TransducerObj.type=='R').*(TransducerObj.platten=='A')));
            n_sT = length(find((TransducerObj.type=='S').*(TransducerObj.platten=='B')));
            n_rT = length(find((TransducerObj.type=='R').*(TransducerObj.platten=='B'))); 
            myMap = zeros((n_sB*n_rB+n_sT*n_rT),2);
            wave_type_mat=zeros((n_sB*n_rB+n_sT*n_rT),2);
    
            pl_f = ['A';'B'];
            k = 1;
            for p = 1:2     
            % find source on platten
                ks= find(((TransducerObj.platten==pl_f(p)) .* (TransducerObj.type=='S'))==1);
                s_type=TransducerObj.wave_mode(ks);
                % find receiver on opposite platten
                kr= find(((TransducerObj.platten==pl_f(p)) .* (TransducerObj.type=='R'))==1) -TransducerObj.n_sources;
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
            
            
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end
        
        
        
        % constructor of diffraction S-R pairs for one side (N-E-S-W)
        % Dong Liu-19/09/2019
        function objpair = SidePairsNeighborPlattens(TransducerObj,platten_list, sidemarker)
            % sidemarker can be 'N', 'S', 'E', 'W' only, representing
            % different sides
            % We only take four transducers on the top and four on the
            % bottom to build the diffracted pairs for the given side. Note 
            % that it is possible to have more transducers on the top and 
            % bottom plattens to build the pairs
            
            n_p=length(platten_list);
            pl_f=char(zeros(n_p,1));
            for p=1:n_p
                pl_f(p)= platten_list(p).face;
            end
            
            if isempty(intersect(pl_f,sidemarker))
                disp(' error in given platten face or corresponding platten list input!');
                return
            end
            
            p_side = find(pl_f==sidemarker) ;
            
            % find source on platten
            ks_1= find(((TransducerObj.platten==platten_list(p_side).id) .* (TransducerObj.type=='S'))==1);
            
            % find receiver on platten
            kr_2= find(((TransducerObj.platten==platten_list(p_side).id) .* (TransducerObj.type=='R'))==1) -TransducerObj.n_sources;
            
            switch sidemarker
                case 'N'
                    knb=1:4;
                    %knb = 9:16;
                    % one can add more transducers here if necessary, the values should be channel values plus one.
                case 'S'
                    knb=5:8;
                    %knb=9:16;
                case 'E'
                    knb=9:12;
                    %knb=1:8;
                case 'W'
                    knb=13:16;
                    %knb=1:8;
                otherwise
                    disp(' error in given platten face or corresponding platten list input!')
            end
            ks_tb=zeros(length(knb),1);
            kr_tb=ks_tb;
            for ki=1:length(knb)
                disp(find(((TransducerObj.channel==knb(ki)-1) .* (TransducerObj.type=='S'))==1));
                ks_tb(ki,1)=find(((TransducerObj.channel==knb(ki)-1) .* (TransducerObj.type=='S'))==1);
                kr_tb(ki,1)=find(((TransducerObj.channel==knb(ki)-1) .* (TransducerObj.type=='R'))==1) -TransducerObj.n_sources;
            end

            s_type=TransducerObj.wave_mode(ks_1);
            r_type=TransducerObj.wave_mode(kr_2);
            stb_type=TransducerObj.wave_mode(ks_tb);
            rtb_type=TransducerObj.wave_mode(kr_tb);
            
            myMap=zeros((length(ks_1) )*(length(kr_tb))+(length(ks_tb) )*(length(kr_2)),2);
            wave_type_mat=zeros((length(ks_1) )*(length(kr_tb))+(length(ks_tb) )*(length(kr_2)),2);
            
            % Source on the side platten and receiver on the top-bottom
            % platten
            k=1;
            for i=1:length(ks_1)
                myMap(k:k+length(kr_tb)-1,1)=ks_1(i);
                myMap(k:k+length(kr_tb)-1,2)=kr_tb;
                wave_type_mat(k:k+length(kr_tb)-1,1)=s_type(i);
                wave_type_mat(k:k+length(kr_tb)-1,2)=rtb_type;
                k=k+length(kr_tb);
            end
            
            % Source on the top-bottom platten and receiver on the side
            % platten
            for j=1:length(ks_tb)
                myMap(k:k+length(kr_2)-1,1)=ks_tb(j);
                myMap(k:k+length(kr_2)-1,2)=kr_2;
                wave_type_mat(k:k+length(kr_2)-1,1)=stb_type(j);
                wave_type_mat(k:k+length(kr_2)-1,2)=r_type;
                k=k+length(kr_2);
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
            
            objpair = SourceReceiverPairs(TransducerObj,platten_list,myMap,wave_type);
            
        end                
        
    end
    
end

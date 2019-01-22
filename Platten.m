classdef Platten
    % A class for the definition of platten and their position during an experiment
    %   Detailed explanation goes here
    
    properties
        
        type ; % cross or grid
        id ;   % id of the block
       
        xy_holes;
        face;
        xloc;  % define the local platten coordinates in the global Block system
        yloc;
        
        offset_x=0.;
        offset_y=0.;
        
    end
    
    
    properties (Dependent)
         n_holes;
    end
    
    
    methods
        
        % constructor
        function obj = Platten(ID,face,xloc,yloc,varargin)
            
            obj.id=ID; % platten id
            obj.face=face; % touching face of the block N,S,E,W,T,B
            obj.xloc=xloc;
            obj.yloc=yloc;
            
            switch ID
                case {'A','B','C'}
                    % it's a cross
                    obj.type='cross';
                    
                    obj.xy_holes=[ 0 0  ; 
                                   0.5 0.5;
                                   1. 1.;
                                   -0.5 0.5;
                                   -1 1;];  % matrix nHoles times 2;

                    
                case {'D'}
                    % it's a grid
                    
                    obj.type='grid';
                    
                    obj.xy_holes=[ 0 0 ];  % matrix nHoles times 2;
                   
                    
            end
            
            % add varargin
            nVarargs = length(varargin);
            if (nVarargs==2)
                obj.offset_x=varargin{1};
                obj.offset_y=varargin{2};
            end
            
            
        end
        % end constructor
        
        function value=get.n_holes(obj)
            value = size(obj.xy_holes,1);
        end
        
        % METHODS
        % 1 method -> plot platten geometry
        
    end
    
end


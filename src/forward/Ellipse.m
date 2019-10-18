% Class for an ellipse
% 
% Brice Lecampion 
classdef Ellipse
    
    properties
        
        % center coordinates
        xc;
        yc;
        zc;
        
        % major and minor axis
        a;
        b;
        
        % euler angles
        Ealpha;
        Ebeta;
        Egamma;
        
        
    end
    
    properties (Dependent)
        
        % Euler rotation matrix zxz
        EulerRot ;
        
    end
    
    
    methods
        
        % constructor
        function obj=Ellipse(a,b, XYZ_c,Ealpha,Ebeta,Egamma)
            
            %  if a>=b
            obj.a=a;
            obj.b=b;
            %             else
            %                 obj.a=b;
            %                 obj.b=a;
            %             end
            
            obj.xc=XYZ_c(1);
            obj.yc=XYZ_c(2);
            obj.zc=XYZ_c(3);
            
            obj.Ealpha=Ealpha;
            obj.Ebeta=Ebeta;
            obj.Egamma=Egamma;
            
        end
        
        % Rotation matrix
        function EulerRot=get.EulerRot(obj)
            
            EulerRot = [cos(obj.Ealpha)*cos(obj.Egamma)-cos(obj.Ebeta)*sin(obj.Ealpha)*sin(obj.Egamma),...
                -cos(obj.Ealpha)*sin(obj.Egamma)-cos(obj.Ebeta)*cos(obj.Egamma)*sin(obj.Ealpha),...
                sin(obj.Ealpha)*sin(obj.Ebeta);...
                cos(obj.Egamma)*sin(obj.Ealpha)+cos(obj.Ealpha)*cos(obj.Ebeta)*sin(obj.Egamma),...
                cos(obj.Ealpha)*cos(obj.Ebeta)*cos(obj.Egamma)-sin(obj.Ealpha)*sin(obj.Egamma),...
                -cos(obj.Ealpha)*sin(obj.Ebeta);...
                sin(obj.Ebeta)*sin(obj.Egamma), cos(obj.Egamma)*sin(obj.Ebeta), cos(obj.Ebeta)];
        end
        
        
        % Get points on Ellipse
        function [Elli_pts]=PointsOnEllipse(obj,npts)
            t=linspace(0.,2*pi,npts+1);
            t=t(1:end-1);
            
            Elli_pts = [obj.a*cos(t);  obj.b*sin(t);  0.*t];
            
            % rotate
            Elli_pts =( [obj.xc+0.*t; obj.yc+0.*t; obj.zc+0.*t]+(obj.EulerRot*Elli_pts) )';
            
        end
        
        % get normal to ellipse
        function [En]=Normal(obj)
            t=[0.,pi/2.];
            Elli_pts =[obj.a*cos(t);  obj.b*sin(t);  0.*t];
            
            Elli_pts =( [obj.xc+0.*t; obj.yc+0.*t; obj.zc+0.*t]+(obj.EulerRot*Elli_pts) )';
            
            
            En = cross(Elli_pts(1,:),Elli_pts(2,:));
            En=En/norm(En);
        end
        
        % plot Ellipse
        function [fig_handle]=plotEllipse(obj,varargin)
            
            narg = length(varargin);
            pl_style='r-';
            if ~isempty(varargin)
                if isgraphics(varargin{1})
                    fig_handle = figure(varargin{1});
                    if narg>1
                        disp(varargin{2});
                        if ischar(varargin{2})
                            pl_style=varargin{2};
                        end
                    end
                    
                else
                    if ischar(varargin{1})
                        pl_style=varargin{1};
                        
                    end
                    fig_handle = figure;
                end
            else
                
                fig_handle = figure;
                
            end
            hold on
            
            %---
            [Elli_pts]=PointsOnEllipse(obj,120); % discretizing the ellipse with 120 pts
            
            plot3(Elli_pts(:,1),Elli_pts(:,2),Elli_pts(:,3),pl_style,'LineWidth',3);
            
            
        end
        
        
        
    end
    
    
end
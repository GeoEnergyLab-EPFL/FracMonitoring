% Class for a radial geometry
% 
% corrected from Ellipse.m
% Dong Liu -- 21/09/2019
classdef Radial
    
    properties
        
        % center coordinates
        xc;
        yc;
        zc;
        
        % radius
        r;

        % euler angles
        Ealpha;
        Ebeta;        
        
    end
    
    properties (Dependent)
        
        % Euler rotation matrix zxz
        EulerRot ;
        
    end
    
    
    methods
        
        % constructor
        function obj=Radial(r, XYZ_c,Ealpha,Ebeta)
            
            obj.r=r;
            
            obj.xc=XYZ_c(1);
            obj.yc=XYZ_c(2);
            obj.zc=XYZ_c(3);
            
            obj.Ealpha=Ealpha;
            obj.Ebeta=Ebeta;
            
        end
        
        % Rotation matrix
        function EulerRot=get.EulerRot(obj)
            
            EulerRot = [cos(obj.Ealpha),...
               -cos(obj.Ebeta)*sin(obj.Ealpha),...
                sin(obj.Ealpha)*sin(obj.Ebeta);...
                sin(obj.Ealpha),...
                cos(obj.Ealpha)*cos(obj.Ebeta),...
                -cos(obj.Ealpha)*sin(obj.Ebeta);...
                0, sin(obj.Ebeta), cos(obj.Ebeta)];
        end
        
        
        % Get points on Radial
        function [Raidal_pts]=PointsOnRadial(obj,npts)
            t=linspace(0.,2*pi,npts+1);
            t=t(1:end-1);
            
            Raidal_pts = [obj.r*cos(t);  obj.r*sin(t);  0.*t];
            
            % rotate
            Raidal_pts =( [obj.xc+0.*t; obj.yc+0.*t; obj.zc+0.*t]+(obj.EulerRot*Raidal_pts) )';
            
        end
        
        % get normal to radial
        function [En]=Normal(obj)
            t=[0.,pi/2.];
            Radial_pts =[obj.r*cos(t);  obj.r*sin(t);  0.*t];
            
            Radial_pts =( [obj.xc+0.*t; obj.yc+0.*t; obj.zc+0.*t]+(obj.EulerRot*Radial_pts) )';
            
            
            En = cross(Radial_pts(1,:),Radial_pts(2,:));
            En=En/norm(En);
        end
        
        % plot Radial
        function [fig_handle]=plotRadial(obj,varargin)
            
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
            [Radial_pts]=PointsOnRadial(obj,120); % discretizing the ellipse with 120 pts
            
            plot3(Radial_pts(:,1),Radial_pts(:,2),Radial_pts(:,3),pl_style,'LineWidth',3);
            
            
        end
        
        
        
    end
    
    
end
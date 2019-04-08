classdef TISolid < Solid
    % description
    % 
    % classs for Transversely isotropic solid properties
    %
    % todo :
    %       Definition of the material frame in the global frame must be
    %       added in the properties (and the dependent TI direction vector
    %       etc.)
    %       Group Velocity - compute Sh & Sv and do some checks
    % 
    properties
        
        % note here all in SI 

        c11;
        c33;
        c13;
        c55;
        c66;
        c44;
        
        % normal to the isotropy plan ein the global system of coordinates 
        
    end
    
    properties (Dependent)
 
        C_ij;%
        
    end
    
    methods
        
        % constructor
        function obj=TISolid(density,c11,c33,c13,c55,c66,varargin)
            
            obj=obj@Solid(density,varargin);
            
            obj.c11=c11;
            obj.c33=c33;
            obj.c13=c13;
            obj.c55=c55;
            obj.c66=c66;
            
            obj.c44=c55;% I don t understand usually one enter c44 not c55, anyway 
            
        end
        
        
        % Group velocity should move here !
        
        function [Vp, Vsv, Vsh, Psi] = GroupVel(obj,theta)
            %GROUPVEL computes the group velocities of a TI medium for discrete angles
            %   GROUPVEL returns the group velocities and group angle of a TI medium as
            %   a function of the elastic constants C11, C33, C13, C55, C66, the
            %   density RHO, for a vector of phase angles theta
            %
            % Thomas Blum, Geo-Energy Lab, EPFL, November 2018
            %
            %  Theta definition ???
            %
            % Equations describing propagation in TI medium from Tsvankin (2001) book
            
            % symbolic derivations of the P-wave group velocity
            % define symbolic variables
            syms thetaSym
            
            tmp1Sym = (obj.c11+obj.c55)*sin(thetaSym).^2 + (obj.c33+obj.c55)*cos(thetaSym).^2;
            tmp2Sym = (obj.c11-obj.c55)*sin(thetaSym).^2 - (obj.c33-obj.c55)*cos(thetaSym).^2;
            tmp3Sym = 4*(obj.c13+obj.c55)^2*sin(thetaSym).^2.*cos(thetaSym).^2;
            
            % compute phase velocity
            % BL: Here I change to SI for consistency & avoiding factors
            VpSym = sqrt((tmp1Sym+sqrt((tmp2Sym).^2+tmp3Sym))/(2*rho));
            
            % and derivative
            dVpSym = diff(VpSym,thetaSym);
            % compute group velocity (Eq 1.70 in Tsvankin)
            UpSym = VpSym.*sqrt(1+(1./VpSym.*dVpSym).^2);
            
            % estimate for a range of phase angles
            Vp = double(subs(UpSym,thetaSym,theta));
            
            % shear velocities (not yet computed)
            Vsh = NaN;
            Vsv = NaN;
            
            % group angle (Eq. 1.71 in Tsvankin)
            num = double(subs(dVpSym/VpSym,thetaSym,theta));
            denom = sin(theta).*cos(theta).*(1-tan(theta).*num);
            tanPsi = tan(theta).*(1+num./denom);
            Psi = atan(tanPsi);
            
            % fix NaN (from division by zero)
            Psi(isnan(Psi)) = 0;
            
            % custom-made unwrap
            dPsi = diff(Psi);
            % fix positive jumps
            Ipos = find(dPsi>=0.9*pi);
            for ii = Ipos
                Psi(ii+1:end) = Psi(ii+1:end)-pi;
            end
            % fix negative jumps
            Ineg = find(dPsi<=-0.9*pi);
            for ii = Ineg
                Psi(ii+1:end) = Psi(ii+1:end)+pi;
            end
            
        end
        
               % get Cij in Voigt notation
        function Cij = get.C_ij(obj)
            
           
            c12=obj.c11-2*obj.c66;
            
            Cij=[obj.c11 c12 obj.c13 0 0 0 ;...
                c12 obj.c11 obj.c13 0 0 0; ...
                obj.c13 obj.c13 obj.c33 0 0 0 ; ...
                0 0 0 obj.c44 0 0 ;...
                0 0 0 0 obj.c44 0 ;...
                0 0 0 0 0 obj.c66];
            
        end
        
    end
end

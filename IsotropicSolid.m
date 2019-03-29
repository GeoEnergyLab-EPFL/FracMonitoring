classdef IsotropicSolid < Solid
% class for an isotropic solid
    properties
        
        Young;
        nu;
        
    end

    properties (Dependent)
        
        K; % bulk modulus
        
        G; % shear modulus

        Vp; % P wave velocity
        
        Vs; % shear wave velocity
        
        C_ij; % Stiffness matrix in Voigt notation
        
    end
    
    methods
    
        % constructor
        function obj=IsotropicSolid(density,E,nu,varargin)
            
            obj=obj@Solid(density,varargin);
            
            obj.Young=E;
            obj.nu=nu;
            
        end
        
         function K=get.K(obj)
            K=obj.Young/(3*(1-2.*obj.nu));
         end
         
         function G=get.G(obj)
            G=obj.Young/(2*(1+obj.nu));
         end
         
         function Vp=get.Vp(obj)
             
            Vp=sqrt((obj.K+2./3.*obj.G)/obj.density);
            
         end
         
         function Vs=get.Vs(obj)
            
              Vs=sqrt(obj.G/obj.density);
            
         end
         
         function Cij = get.C_ij(obj)
         
            La=obj.K+(4./3.)*obj.G;
            Lb=obj.K-(2./3.)*obj.G;

            Cij=[La Lb Lb 0 0 0 ;...
                Lb La Lb 0 0 0; ...
                Lb Lb La 0 0 0 ; ...
                0 0 0 2*obj.G 0 0 ;...
                0 0 0 0 2*obj.G 0 ;...
                0 0 0 0 0 2*obj.G];
                
         end
         
    end
    
    
end

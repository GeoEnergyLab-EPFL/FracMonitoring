classdef GaussianPrior
% gaussian prior pdf class
% only uncorrelated prior - i.e. covariance matrix is diagonal
    properties
    
        mp;
        
        invCpdiag;
        
    end
    
    methods
        
        function [obj]=GaussianPrior(mp,sigmas)
           
            if (length(mp)~= length(sigmas) )
                disp('error - prior values and  variance not of the same length');
                return
            end
            
            obj.mp=mp;
            
            obj.invCpdiag = 1./(sigmas.^2) ; % diagonal matrix - we only store the diagonal

        end
        
        function [mPriorPdf]=minusLogPrior(obj,z)
            %  -log Prior 
            % 
            % Gaussian prior
            %   
            %   
            m=z(1:length(z)-1);
            
            mPriorPdf=0.5*((m-obj.mp)'*((m-obj.mp).*(obj.invCpdiag)));
            
        end
        
   end
   
end
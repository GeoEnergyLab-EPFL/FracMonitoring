function [mpost,Ctilde]=PostProcessInversion(inv_opt,accept,PX,X,k_o,stab,mcmc_Diag_res_logfile)
%
% INPUT PARAMETER
%  inv_opt :: inversion options structure
%  accept  :: rate of MCMC acceptation
%  X       :: MCMC vector outputs   a matrix number_of_MCMC_steps by number
%               of parameters
%  PX      :: a vector (length number_of_MCMC_steps) containing value
%             of -log Post
%  k_o     :: integer (define the end of the burn in period)
%  stab    :: integer flag (0, or 1) stabilized chain (1) or not (0)
%  mcmc_Diag_res_logfile 
%          :: file name for MCMC diagnostic output
%
% OUTPUT 
%   mpost  :: most probable value of the model parameter vector  ( q by p   MATRIX In
%   the case of p modes )
%   Ctilde :: posterior covariance matrix around mpost. :: q by q by p
%   matrices in the case of p modes.
%

global G measured Cp prior_information;
  
   if (inv_opt.MCMC==1) && (length(PX)>1)
       % MCMC diagnostic
       
       [mpost,Ctilde]=ProcessingMCMC(accept,PX,X,k_o,stab,mcmc_Diag_res_logfile);
       % takes the quadratic approximation covariance matrix for Ctilde for now (ideally user
       % should be alllowed to choose
  
        
   else
       % NO MCMC has been performed
       mpost=X;
       [Ctilde]=Posterior_Covariance_Matrix(mpost, G, measured, Cp, prior_information);       
       
   end

return

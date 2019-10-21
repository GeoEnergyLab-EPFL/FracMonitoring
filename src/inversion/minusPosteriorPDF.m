function [mpost]=minusPosteriorPDF(m)
 % wrapper to compute -Log Posterior PDF
 % note that d prior  must be declared as global variable
 % as well as Solid (solid properties object)  SRPairs (source receivers pairs)
 % and diagonal of inverse of measurmenets covariance matrix (Cdinvdiag)
 %
 %
 global prior ;
    % both for the the radial and ellipse case
    mL=minusLikelihoodEllipseDiff(m); % compute -Log Likelihood
    
    mP=minusLogPrior(prior,m);% compute -Log Prior pdf
    
    mpost=mL+mP;
    
end

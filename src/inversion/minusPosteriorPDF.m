function [zpost]=minusPosteriorPDF(z)
 % wrapper to compute -Log Posterior PDF
 % note that d prior  must be declared as global variable
 % as well as Solid (solid properties object)  SRPairs (source receivers pairs)
 % and diagonal of inverse of measurmenets covariance matrix (Cdinvdiag)
 %
 %
 global prior ;
 global d;
    % both for the the radial and ellipse case
    
    % the model parameters only cover part of the vector m
    zL=minusLikelihoodEllipseDiff(z); % compute -Log Likelihood
    
    zP=minusLogPrior(prior,z);% compute -Log Prior pdf
    
    znoise=z(length(z)).*length(d);
    
    zpost=zL+zP+znoise;
    
end

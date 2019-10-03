function [mLikelihood]=minusLikelihoodEllipseDiff(m)

 % wrapper to the forward diffraction problem
 % Gaussian likelihood with uncorrelated measurement covariance
 % 
 % Solid (solid properties object)  SRPairs (source receivers pairs)
 % and diagonal of inverse of measurmenets covariance matrix (Cdinvdiag)
 %
    global Solid SRPairs Cdinvdiag d;
    
    % create ellipse object
    
    % case 8
    if (length(m)~=8)
        disp(' error in number of parameters - ');
        return
    end
    
    ell=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
    
    res= diffractionForward(Solid,SRPairs,ell);

    G=res(:,1);
    
    if (length(G)~=length(d) )
        disp('error in number of data / predictions');
    end
        
    mLikelihood=0.5*((G-d)'*((G-d).*(Cdinvdiag)));
    
end

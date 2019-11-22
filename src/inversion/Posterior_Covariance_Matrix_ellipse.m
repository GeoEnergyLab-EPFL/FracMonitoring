function [Ctilde]=Posterior_Covariance_Matrix_ellipse(mpost)
%
% This function computes a quadratic approximation  of the covariance
% matrix of the posterior pdf 
%  Ctilde=(X'CD^{-1} X +  Cp^{-1} ) ^{-1}
%  The jacobian matrices are computed via finite difference.
%  measurements part (splitted in magnitude and orientation)
%  CHANGE & PUT FLAGS ON MEASURED VALUE
%
% Inputs:
%  mpost             :: 1-dimensional array which contains the inverted 
%                       posterior m
% Outputs:
%  Ctilde            :: (m_i, m_j), 2D array which defines the posterior 
%                       covariance matrix
%

global Solid SRPairs Cdinvdiag d prior ray_type;

global z_c; % for the case of the fixed z-coordinate;


inv_Cd=diag(Cdinvdiag);
inv_Cp=diag(prior.invCpdiag);
  
%%%% FINITE DIFFERENCE APPROXIMATION OF JACOBIAN Matrices

q=length(mpost);

h=1e-4;


X=zeros(length(d),q);

for i=1:q
    
    m_plus=mpost;
    m_minus=mpost;
    m_plus(i)=mpost(i)*(1+h);
    m_minus(i)=mpost(i)*(1-h);

    m=m_plus;
    % case of flexible z_c coordiante
%     if length(m)>6
%         shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
%     else
%         shape=Radial(m(1),m(2:4),m(5),m(6));
%     end

    % case of fixed z_c coordiante
    switch length(m)
        case 8
            shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 7
            shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7));
        case 6 
            shape=Radial(m(1),m(2:4),m(5),m(6));
        case 5
            shape=Radial(m(1),[m(2:3); z_c] ,m(4),m(5));
        otherwise
            disp('Please check your input');
    end
    
    res= diffractionForward(Solid,SRPairs,shape,ray_type);
    G_plus=res(:,1);
    m=m_minus;
    % case of flexible z-coordinate
%     if length(m)>6
%         shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
%     else
%         shape=Radial(m(1),m(2:4),m(5),m(6));
%     end
    % case of fixed z-coordinate 
    switch length(m)
        case 8
            shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 7
            shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7));
        case 6 
            shape=Radial(m(1),m(2:4),m(5),m(6));
        case 5
            shape=Radial(m(1),[m(2:3); z_c] ,m(4),m(5));
        otherwise
            disp('Please check your input');
    end
    
    res= diffractionForward(Solid,SRPairs,shape,ray_type);
    G_minus=res(:,1);
    
    X(:,i)=(G_plus-G_minus)/(2*h*mpost(i));
    
end

Ctilde_minus= X'*inv_Cd*X+inv_Cp;

Ctilde=pinv(Ctilde_minus);


return

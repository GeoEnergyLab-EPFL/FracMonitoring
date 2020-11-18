function [Ctilde]=Posterior_Covariance_Matrix_ellipse(zpost)
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

global m_ind;

global Solid SRPairs d prior ray_type;

global z_c; % for the case of the fixed z-coordinate;

mpost=zpost(1:length(zpost)-1);
noise=zpost(length(zpost));

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if length(m)>6
%         shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
%     else
%         shape=Radial(m(1),m(2:4),m(5),m(6));
%     end

%     % case of fixed z_c coordiante
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     switch length(m)
%         case 8
%             shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
%         case 7
%             shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7));
%         case 6 
%             shape=Radial(m(1),m(2:4),m(5),m(6));
%         case 5
%             shape=Radial(m(1),[m(2:3); z_c] ,m(4),m(5));
%         otherwise
%             disp('Please check your input');
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch m_ind
        case 1
            shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 2
            shape=Radial(m(1),m(2:4),m(5),m(6));
        case 3
            shape=Ellipse(m(1),m(2),m(3:5),m(6),0,0);
        case 4
            shape=Radial(m(1),m(2:4),0,0);
        case 5
            shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7));
        case 6
            shape=Radial(m(1),[m(2:3); z_c] ,m(4),m(5));
        otherwise 
            disp('Please check your input');
    end
    
    res= diffractionForward(Solid,SRPairs,shape,ray_type);
    G_plus=res(:,1);
    m=m_minus;
    % case of flexible z-coordinate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if length(m)>6
%         shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
%     else
%         shape=Radial(m(1),m(2:4),m(5),m(6));
%     end
    % case of fixed z-coordinate 
%     switch length(m)
%         case 8
%             shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
%         case 7
%             shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7));
%         case 6 
%             shape=Radial(m(1),m(2:4),m(5),m(6));
%         case 5
%             shape=Radial(m(1),[m(2:3); z_c] ,m(4),m(5));
%         otherwise
%             disp('Please check your input');
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch m_ind
        case 1
            shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 2
            shape=Radial(m(1),m(2:4),m(5),m(6));
        case 3
            shape=Ellipse(m(1),m(2),m(3:5),m(6),0,0);
        case 4
            shape=Radial(m(1),m(2:4),0,0);
        case 5
            shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7));
        case 6
            shape=Radial(m(1),[m(2:3); z_c] ,m(4),m(5));
        otherwise 
            disp('Please check your input');
    end
    
    
    res= diffractionForward(Solid,SRPairs,shape,ray_type);
    G_minus=res(:,1);
    
    X(:,i)=(G_plus-G_minus)/(2*h*mpost(i));
    
    % calculate G using the mpost
    m=mpost;
    switch m_ind
        case 1
            shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 2
            shape=Radial(m(1),m(2:4),m(5),m(6));
        case 3
            shape=Ellipse(m(1),m(2),m(3:5),m(6),0,0);
        case 4
            shape=Radial(m(1),m(2:4),0,0);
        case 5
            shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7));
        case 6
            shape=Radial(m(1),[m(2:3); z_c] ,m(4),m(5));
        otherwise 
            disp('Please check your input');
    end
    
    
    res= diffractionForward(Solid,SRPairs,shape,ray_type);
    G_mpost=res(:,1);

    

    
end

N_p=length(mpost);
Ctilde_minus=ones(N_p+1);
H_ij=(exp(-2*noise))*(X'*X)+inv_Cp;
H_if=-2*exp(-2*noise)*(X'*(G_mpost-d));
H_ff=2*exp(-2*noise)*((d-G_mpost)'*(d-G_mpost));
for i=1:N_p+1
    if(i<N_p+1)
        for j=1:N_p
            Ctilde_minus(i,j)=H_ij(i,j);
        end
    else
        for j=1:N_p
            Ctilde_minus(N_p+1,j)=H_if(j);
            Ctilde_minus(j,N_p+1)=H_if(j);
        end
        Ctilde_minus(N_p+1,N_p+1)=H_ff;
    end
end

Ctilde=pinv(Ctilde_minus);


return

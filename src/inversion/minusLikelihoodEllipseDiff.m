function [mLikelihood]=minusLikelihoodEllipseDiff(m)

 % wrapper to the forward diffraction problem
 % Gaussian likelihood with uncorrelated measurement covariance
 % 
 % Solid (solid properties object)  SRPairs (source receivers pairs)
 % and diagonal of inverse of measurmenets covariance matrix (Cdinvdiag)
 %
    global Solid SRPairs Cdinvdiag d ray_type;
    
    % for the case of the fixed z_c coordinate
    global z_c;
    
    % create ellipse object
    if(length(m)>6)
%     % case 8
%         if (length(m)~=8)
%             disp(' error in number of parameters - ');
%             return
%         end
%     
%         shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
    % case 7 and case 8 18/11/2019
        if (length(m)>8)
            disp(' error in number of parameters - ');
            return
        end
        if (length(m)==8)
            shape=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        elseif (length(m)==7)
            shape=Ellipse(m(1),m(2),[m(3:4); z_c],m(5),m(6),m(7)); 
            % z_c is the fixed z-cooridinate (~fracture center coordinate)
        end
    
    else
        % create radial object
%         % case 6
%         if (length(m)~=6)
%             disp(' error in number of parameters - ');
%             return
%         end
%         shape=Radial(m(1),m(2:4),m(5),m(6));
        % case 6 and 5 
        if (length(m)<5)
            disp(' error in number of parameters - ');
            return
        end
        if (length(m)==6)
            shape=Radial(m(1),m(2:4),m(5),m(6));
        elseif (length(m)==5)
            shape=Radial(m(1),[m(2:3); z_c],m(4),m(5));
        end
        
    end
    
    
    res= diffractionForward(Solid,SRPairs,shape,ray_type);

    G=res(:,1);
    
    if (length(G)~=length(d) )
        disp('error in number of data / predictions');
    end
        
    mLikelihood=0.5*((G-d)'*((G-d).*(Cdinvdiag)));
    
end

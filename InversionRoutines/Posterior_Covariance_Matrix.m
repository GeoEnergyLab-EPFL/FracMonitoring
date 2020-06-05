function [Ctilde]=Posterior_Covariance_Matrix(mpost, G, measured, Cp, prior_information)
%
% This function computes a quadratic approximation  of the covariance
% matrix of the posterior pdf 
%  Ctilde=(X'CD^{-1} X + Sum Y' Cp^{-1} Y ) ^{-1}
%  The jacobian matrices are computed via finite difference.
%  measurements part (splitted in magnitude and orientation)
%  CHANGE & PUT FLAGS ON MEASURED VALUE
%
% Inputs:
%  mpost             :: 1-dimensional array which contains the inverted 
%                       posterior m
%  G(loading_id, measure_point_id, stress_id) 
%                    :: a 3-dimensional array which contains the stress 
%                   magnitudes at the measurement stations of each 
%                   fundamental loading(stress_id is from 1 to 6 
%                   indicating "XX", "YY", "ZZ", "XY", "YZ", "ZX")
%  measured          :: A structure to define the measured data and 
%                       uncertainties
%  Cp                :: (m_i, m_j), 2D array which defines the prior 
%                       covariance matrix
%  prior_information :: A structure to define the prior information data 
%                       and uncertainties
% Outputs:
%  Ctilde            :: (m_i, m_j), 2D array which defines the posterior 
%                       covariance matrix
%

Ns_prior=prior_information.Ns;

Ns=measured.Ns;


%wflag_mag=find(measured.mag_flag~=0);
wflag_or=find(measured.orient_flag~=0);

%data_mag=measured.mag_data(wflag_mag);
data_orient=measured.orient_data(wflag_or,:);
num_data_orient=length(data_orient(:,1));

%%%% covariance matrix on stresses magnitudes.
%Cd_mag=diag(measured.mag_uncertainty(wflag_mag).^2);
%Cd_orient=diag((measured.orient_uncertainty(wflag_or)).^2);

m=mpost;

[d_magnitude,d_orientation]=Forward_Stress(m,G,Ns,measured.mag_flag,measured.orient_flag);
%%compute orientation function:: dot products 
if(num_data_orient>0)
    v=zeros(length(data_orient(:,1)),1);
    for i=1:length(data_orient(:,1)),
        index=(i-1)+1;
        v(i)=180.0*my_acos(abs(d_orientation(index,:)*(data_orient(index,:)')))/pi;
    end;
end

[inv_Cd_mag,inv_Cd_orient,inv_Cp_sv, inv_Cp_sh,inv_Cp_shmin_sv_ratio,inv_Cp_shmax_sv_ratio]=InitializeCova(measured, prior_information, Ns_prior);

%%%%%%%% GEOLOGICAL PRIORS
[I_SV, I_SH,I_shmin_sv_ratio, I_shmax_sv_ratio]=IGeological_Prior(m);


%%%% FINITE DIFFERENCE APPROXIMATION OF JACOBIAN Matrices

q=length(mpost);

h=1e-6;

Xmag=zeros(length(d_magnitude),q);
if(num_data_orient>0)
    Xorient=zeros(length(v),q);
end
%Yregime=zeros(length(I_regime),q);
YSVorient=zeros(length(I_SV),q);
YSHorient=zeros(length(I_SH),q);
Yshmin_sv=zeros(length(I_shmin_sv_ratio),q);
Yshmax_sv=zeros(length(I_shmax_sv_ratio),q);

for i=1:q,
    m_plus=mpost;
    m_minus=mpost;
    m_plus(i)=mpost(i)*(1+h);
    dx_plus=h*mpost(i);
    m_minus(i)=mpost(i)*(1-h);
    dx_minus=h*mpost(i);
    
    [dplus_magnitude,dplus_orientation]=Forward_Stress(m_plus,G,Ns,measured.mag_flag,measured.orient_flag);
    [dminus_magnitude,dminus_orientation]=Forward_Stress(m_minus,G,Ns,measured.mag_flag,measured.orient_flag);
    
    %%compute orientation function:: dot products 
    if(num_data_orient>0)
        v_plus=zeros(length(data_orient(:,1)),1);
        v_minus=zeros(length(data_orient(:,1)),1);
        for j=1:length(data_orient(:,1)),
            index=(j-1)+1;
            v_plus(j)=180.0*my_acos(abs(dplus_orientation(index,:)*(data_orient(index,:)')))/pi;
            v_minus(j)=180.0*my_acos(abs(dminus_orientation(index,:)*(data_orient(index,:)')))/pi;
        end;
        Xorient(:,i)=(v_plus-v_minus)/(dx_plus+dx_minus);
    end
    Xmag(:,i)=(dplus_magnitude-dminus_magnitude)/(dx_plus+dx_minus);

    
    % Geological priors
    [Iplus_SV, Iplus_SH,Iplus_shmin_sv_ratio, Iplus_shmax_sv_ratio]=IGeological_Prior(m_plus);
    [Iminus_SV, Iminus_SH,Iminus_shmin_sv_ratio, Iminus_shmax_sv_ratio]=IGeological_Prior(m_minus);
    %Yregime(:,i)=(Ih_regime-I_regime)/dx;
    if(length(I_SV)>0)
      YSVorient(:,i)=(Iplus_SV-Iminus_SV)/(dx_plus+dx_minus);
    end
    if(length(I_SH)>0)
      YSHorient(:,i)=(Iplus_SH-Iminus_SH)/(dx_plus+dx_minus);
    end
    if (length(I_shmin_sv_ratio)>0);
      Yshmin_sv(:,i)=(Iplus_shmin_sv_ratio-Iminus_shmin_sv_ratio)/(dx_plus+dx_minus);    
    end;
    if (length(I_shmax_sv_ratio)>0);
      Yshmax_sv(:,i)=(Iplus_shmax_sv_ratio-Iminus_shmax_sv_ratio)/(dx_plus+dx_minus);    
    end;
end,

Ctilde_minus= 0.0 ... 
    +Xmag'*inv_Cd_mag*Xmag ...
    +YSVorient'*inv_Cp_sv*YSVorient...
    +YSHorient'*inv_Cp_sh*YSHorient...
    +Yshmin_sv'*inv_Cp_shmin_sv_ratio*Yshmin_sv ...
    +Yshmax_sv'*inv_Cp_shmax_sv_ratio*Yshmax_sv...
    +inv(Cp) ;
if(num_data_orient>0)
    Ctilde_minus= Ctilde_minus +Xorient'*inv_Cd_orient*Xorient;
end
Ctilde=pinv(Ctilde_minus);


return

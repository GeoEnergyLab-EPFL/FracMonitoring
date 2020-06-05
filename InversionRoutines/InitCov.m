function InitCov(mprior)

%
% function to initialize covariance matrix Cp, inv_Cd_* (2 items), inv_Cp_* (5 items)
%
% Input
%	mprior		prior m (vector)
%
% Output
%	none		(saved as global variables)
% 

    global prior_information measured ;

    global Cp ;
    global inv_Cd_mag inv_Cd_orient;
    global inv_Cp_sv inv_Cp_sh inv_Cp_shmin_sv_ratio inv_Cp_shmax_sv_ratio;
    global inv_Cp_Q;
    

    q=length(mprior);
    stprior=abs(mprior)*100; % "flat" like prior
    Cp=diag(stprior.^2);

%init  inv_Cd_mag,inv_Cd_orient,inv_Cp_sv,inv_Cp_sh,inv_Cp_shmin_sv_ratio,inv_Cp_shmax_sv_ratio (set as globals)
    %[inv_Cd_mag,inv_Cd_orient,inv_Cp_sv,inv_Cp_sh,inv_Cp_shmin_sv_ratio,inv_Cp_shmax_sv_ratio]=InitializeCova(measured, prior_information, prior_information.Ns);
    [inv_Cd_mag,inv_Cd_orient,inv_Cp_sv,inv_Cp_sh,inv_Cp_shmin_sv_ratio,inv_Cp_shmax_sv_ratio,inv_Cp_Q]=InitializeCova(measured, prior_information, prior_information.Ns);

return;

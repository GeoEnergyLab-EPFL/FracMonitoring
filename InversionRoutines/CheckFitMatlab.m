function []=CheckFitMatlab(mym,Ccov)
%
% PURPOSE : 
%           plot the fit for each magnitude, separating each well
%
% INPUT PARAMETERS
%   mym  :: the value of the param vector 
%   Ccov :: the corresponding posterior covariance matrix
%

 global G  measured
  
 
 stress_factor=1e-6;
 
 Ns=size(measured.location,1);
  
 [magnitude_fwd,V_magnitude_fwd,orientation_fwd,U_orientation_fwd]=Forward_Stress_and_Uncertainty(mym,Ccov,G, Ns); %currently mag_flag and orient_flag are not used
    
 %%%   wells
 totwellid=[measured.id measured.id measured.id]';
 totwelldepth=-[measured.location(:,3) measured.location(:,3) measured.location(:,3) ]';
 % flag for measured magnitude.
 
    
  %---- check the different components...
  sv_index=[1:3:length(measured.mag_flag)]';
  sh_index=[2:3:length(measured.mag_flag)]';
  sH_index=[3:3:length(measured.mag_flag)]';  
  
  sv_mag_flag=find(measured.mag_flag(1:3:length(measured.mag_flag))==1);
  sh_mag_flag=find(measured.mag_flag(2:3:length(measured.mag_flag))==1);
  sH_mag_flag=find(measured.mag_flag(3:3:length(measured.mag_flag))==1);
 
  %--- measured.mag_data index for different components 
  k_sH_mag=sH_index(sH_mag_flag);
  k_sh_mag=sh_index(sh_mag_flag);
  k_sv_mag=sv_index(sv_mag_flag);
   
  sv_meas_mag=measured.mag_data(k_sv_mag); w_Id_sv=totwellid(k_sv_mag); depth_sv=totwelldepth(k_sv_mag);
  sh_meas_mag=measured.mag_data(k_sh_mag); w_Id_sh=totwellid(k_sh_mag); depth_sh=totwelldepth(k_sh_mag);
  sH_meas_mag=measured.mag_data(k_sH_mag); w_Id_sH=totwellid(k_sH_mag); depth_sH=totwelldepth(k_sH_mag);
    
  %--- corresponding... forward results d_magnitude...
  
  sv_fwd_mag=magnitude_fwd(k_sv_mag); 
  sh_fwd_mag=magnitude_fwd(k_sh_mag);
  sH_fwd_mag=magnitude_fwd(k_sH_mag);
  
  % and std deivation
  Usv_fwd_mag=V_magnitude_fwd(k_sv_mag);
  Ush_fwd_mag=V_magnitude_fwd(k_sh_mag);
  UsH_fwd_mag=V_magnitude_fwd(k_sH_mag);

  
%%%%%%%%%%%%%%%%%%%%%%%%%% ORIENTATION
for i=1:Ns
      
      j=(i-1)*3+1;
 
      SvShminShmax_Vect=orientation_fwd(j:j+2,1:3)';
      
      svaz_hat(i)=atan2(SvShminShmax_Vect(1,1), SvShminShmax_Vect(2,1))/pi*180; %azimuth from North(Y axis)
      shminaz_hat(i)=atan2(SvShminShmax_Vect(1,2), SvShminShmax_Vect(2,2))/pi*180; %azimuth from North(Y axis)
      shmaxaz_hat(i)=atan2(SvShminShmax_Vect(1,3), SvShminShmax_Vect(2,3))/pi*180; %azimuth from North(Y axis)

      if (svaz_hat(i)<0)
        svaz_hat(i)=svaz_hat(i)+180;
      end
      
      if  ( shminaz_hat(i)  <0)
           shminaz_hat(i)=shminaz_hat(i)+180;
      end
      
      if  ( shmaxaz_hat(i)  <0)
           shmaxaz_hat(i)=shmaxaz_hat(i)+180;
      end
      
      
      svpl_hat(i)=abs(atan2(SvShminShmax_Vect(3,1), norm([SvShminShmax_Vect(1,1),SvShminShmax_Vect(2,1)]))/pi*180.0); %atan(z,norm(x,y))
      shminpl_hat(i)=abs(atan2(SvShminShmax_Vect(3,2), norm([SvShminShmax_Vect(1,2),SvShminShmax_Vect(2,2)]))/pi*180.0);
      shmaxpl_hat(i)=abs(atan2(SvShminShmax_Vect(3,3), norm([SvShminShmax_Vect(1,3),SvShminShmax_Vect(2,3)]))/pi*180.0);

   
      %%% U 
       SvShminShmax_Vect=U_orientation_fwd(j:j+2,1:3)';
      
      U_svaz_hat=atan2(SvShminShmax_Vect(1,1), SvShminShmax_Vect(2,1))/pi*180; %azimuth from North(Y axis)
      U_shminaz_hat=atan2(SvShminShmax_Vect(1,2), SvShminShmax_Vect(2,2))/pi*180; %azimuth from North(Y axis)
      U_shmaxaz_hat=atan2(SvShminShmax_Vect(1,3), SvShminShmax_Vect(2,3))/pi*180; %azimuth from North(Y axis)

      if  (U_svaz_hat  <0)
          U_svaz_hat=U_svaz_hat+180;
      end

       if  (  U_shminaz_hat <0)
           U_shminaz_hat=U_shminaz_hat+180;
       end

       if  (  U_shmaxaz_hat <0)
           U_shmaxaz_hat=U_shmaxaz_hat+180;
       end
     
      U_svpl_hat=abs(atan2(SvShminShmax_Vect(3,1), norm([SvShminShmax_Vect(1,1),SvShminShmax_Vect(2,1)]))/pi*180.0); %atan(z,norm(x,y))
      U_shminpl_hat=abs(atan2(SvShminShmax_Vect(3,2), norm([SvShminShmax_Vect(1,2),SvShminShmax_Vect(2,2)]))/pi*180.0);
      U_shmaxpl_hat=abs(atan2(SvShminShmax_Vect(3,3), norm([SvShminShmax_Vect(1,3),SvShminShmax_Vect(2,3)]))/pi*180.0);
                  
      %%%% STD DEVIATION
      std_svaz(i)=abs(svaz_hat(i)-U_svaz_hat);
      std_shminaz(i)=abs(shminaz_hat(i)-U_shminaz_hat);
      std_shmaxaz(i)=abs(shmaxaz_hat(i)-U_shmaxaz_hat);
      
      std_svpl(i)=abs(svpl_hat(i)-U_svpl_hat);
      std_shminpl(i)=abs(shminpl_hat(i)-U_shminpl_hat);
      std_shmaxpl(i)=abs(shmaxpl_hat(i)-U_shmaxpl_hat);
           
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% switch back measurement into plunge etc....
      
      SvShminShmax_Vect=measured.orient_data(j:j+2,1:3)';
      
      svaz_mea(i)=atan2(SvShminShmax_Vect(1,1), SvShminShmax_Vect(2,1))/pi*180; %azimuth from North(Y axis)
      
      shminaz_mea(i)=atan2(SvShminShmax_Vect(1,2), SvShminShmax_Vect(2,2))/pi*180; %azimuth from North(Y axis)
      shmaxaz_mea(i)=atan2(SvShminShmax_Vect(1,3), SvShminShmax_Vect(2,3))/pi*180; %azimuth from North(Y axis)

       if (svaz_mea(i)<0)
            svaz_mea(i)=svaz_mea(i)+180;
       end
      
       if  ( shminaz_mea(i)  <0)
           shminaz_mea(i)=shminaz_mea(i)+180;
       end
      
       if  ( shmaxaz_hat(i)  <0)
           shmaxaz_mea(i)=shmaxaz_mea(i)+180;
       end
    
      svpl_mea(i)=abs(atan2(SvShminShmax_Vect(3,1), norm([SvShminShmax_Vect(1,1),SvShminShmax_Vect(2,1)]))/pi*180.0); %atan(z,norm(x,y))
      shminpl_mea(i)=abs(atan2(SvShminShmax_Vect(3,2), norm([SvShminShmax_Vect(1,2),SvShminShmax_Vect(2,2)]))/pi*180.0);
      shmaxpl_mea(i)=abs(atan2(SvShminShmax_Vect(3,3), norm([SvShminShmax_Vect(1,3),SvShminShmax_Vect(2,3)]))/pi*180.0);
       
     
       %%% MEAS VARIANCE.....
      %std_shaz_m(i)=measured.orient_uncertainty(j+2); % may be should figure out the uncertainty diff ....
      %std_sHaz_m(i)=measured.orient_uncertainty(j+1); 
      %std_svaz_m(i)=measured.orient_uncertainty(j); 
      std_svaz_m(i)=measured.orient_uncertainty(j); 
      std_shaz_m(i)=measured.orient_uncertainty(j+1);
      std_sHaz_m(i)=measured.orient_uncertainty(j+2); 
      
      %std_svpl_m(i)=measured.orient_uncertainty(j+1); 
      %std_shminpl_m(i)=measured.orient_uncertainty(j+1); 
      %std_shmaxpl_m(i)=measured.orient_uncertainty(j+1); 
      std_svpl_m(i)=measured.orient_uncertainty(j); 
      std_shminpl_m(i)=measured.orient_uncertainty(j+1); 
      std_shmaxpl_m(i)=measured.orient_uncertainty(j+2); 
      
      
      % measured.orient_flag(j) Sv (plunge and azimuth)       
      % measured.orient_flag(j+1) Sh      
      % measured.orient_flag(j+2) SH
      
      if ((measured.orient_flag(j)==0) ||(measured.orient_flag(j+1)==0))
        flag_sh(i)=0;
      else
         flag_sh(i)=1;         
      end
      
      if ((measured.orient_flag(j)==0) ||(measured.orient_flag(j+2)==0))
        flag_sH(i)=0;
      else
         flag_sH(i)=1;         
      end
      
      if ((measured.orient_flag(j)==0))
          flag_sv(i)=0;
          
      else
          flag_sv(i)=1;
      end
      
      
end
  
  ksv_O=find(flag_sv==1);
  ksh_O=find(flag_sh==1);
  ksH_O=find(flag_sH==1); 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                 PLOTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% NUMBER OF WELLS
  n_wells=length(unique(measured.id));

  % look for minimum and maximum depth.
  d_max=max(measured.location(:,3));
  d_min=min(measured.location(:,3));
  
  
  %%% SV PLOT for each welll
  handler_sv=figure('Name','SV for each well' );
  d_max=max(depth_sv);
  d_min=min(depth_sv);

  for i=1:n_wells,
     
       jaux=find(w_Id_sv(:)==i);
      
      subplot(1,n_wells,i)
       
       if (~isempty(jaux))
       herrorbar(stress_factor*sv_meas_mag(jaux),depth_sv(jaux),stress_factor*measured.mag_uncertainty(k_sv_mag(jaux)),'.r'); hold on;  
       herrorbar(stress_factor*sv_fwd_mag(jaux),depth_sv(jaux),stress_factor*Usv_fwd_mag(jaux),'.b'); 
       ylim([d_min d_max])
       set(gca,'yDir','reverse');
       end
       title(['Well ID# ', num2str(i)]);
       
  end
  if ( ~isempty(sv_meas_mag) && ~isempty(jaux) )
    %legend(' ','Data',' ', 'Fit ', 'location','Best')
    legend(' ','Measurement Data',' ', 'Inversion Results', 'location','Best')
  else
    fprintf(1, '\n');
    fprintf(2, 'Warning: \n');
    fprintf(1, '	no measurement data for Sv\n');
    fprintf(1, '	empty plot will be generated\n');
    %fprintf(1, '\n');
  end

if ( ~isempty(sv_meas_mag) )
  fprintf('\n');
  fprintf('Sv\n');
  fprintf('mea_mag  mea_std    fit_mag  fit_std\n');
  for i=1:n_wells,
       jaux=find(w_Id_sv(:)==i);
       
       if (~isempty(jaux))
	 mea_mag = stress_factor*sv_meas_mag(jaux);
	 mea_std = stress_factor*measured.mag_uncertainty(k_sv_mag(jaux));
	 fit_mag = stress_factor*sv_fwd_mag(jaux);
	 fit_std = stress_factor*Usv_fwd_mag(jaux);
	 fprintf('%8.2f %8.2f %8.2f %8.2f\n', [mea_mag(:), mea_std(:), fit_mag(:), fit_std(:)]');
       end
  end
end


  %%% SH 
    handler_shmax=figure('Name','SH for each well' );
  d_max=max(depth_sH);
  d_min=min(depth_sH);

  for i=1:n_wells,
        jaux=find(w_Id_sH(:)==i);
      
      subplot(1,n_wells,i)
        if (~isempty(jaux))    
        herrorbar(stress_factor*sH_meas_mag(jaux),depth_sH(jaux),stress_factor*measured.mag_uncertainty(k_sH_mag(jaux)),'.r'); hold on;  
        herrorbar(stress_factor*sH_fwd_mag(jaux),depth_sH(jaux),stress_factor*UsH_fwd_mag(jaux),'.b');
       ylim([d_min d_max])
       set(gca,'yDir','reverse');

        end
       title(['Well ID# ', num2str(i)]);
  end
  if ( ~isempty(sH_meas_mag) && ~isempty(jaux) )
    %legend(' ','Data',' ', 'Fit ', 'location','Best')
    legend(' ','Measurement Data',' ', 'Inversion Results', 'location','Best')
  else
    fprintf(1, '\n');
    fprintf(2, 'Warning: \n');
    fprintf(1, '	no measurement data for shmax\n');
    fprintf(1, '	empty plot will be generated\n');
    %fprintf(1, '\n');
  end
 

  % Sh, 3 wells raw plots
    handler_shmin=figure('Name','Sh for each well' );
      d_max=max(depth_sh);
  d_min=min(depth_sh);

  for i=1:n_wells,
      jaux=find(w_Id_sh(:)==i);
       
      subplot(1,n_wells,i)
      if (~isempty(jaux))         
       herrorbar(stress_factor*sh_meas_mag(jaux),depth_sh(jaux),stress_factor*measured.mag_uncertainty(k_sh_mag(jaux)),'.r'); hold on;  
       herrorbar(stress_factor*sh_fwd_mag(jaux),depth_sh(jaux),stress_factor*Ush_fwd_mag(jaux),'.b'); 
       ylim([d_min d_max])
       xlim([0 40])
       set(gca,'yDir','reverse');
      end 
       title(['Well ID# ', num2str(i)]);
  end
  if ( ~isempty(sh_meas_mag) && ~isempty(jaux) )
    %legend(' ','Data',' ', 'Fit ', 'location','Best')
    legend(' ','Measurement Data',' ', 'Inversion Results', 'location','Best')
  else
    fprintf(1, '\n');
    fprintf(2, 'Warning: \n');
    fprintf(1, '	no measurement data for shmin\n');
    fprintf(1, '	empty plot will be generated\n');
    %fprintf(1, '\n');
  end
if ( ~isempty(sh_meas_mag) )
  fprintf('\n');
  fprintf('Shmin\n');
  fprintf('mea_mag  mea_std    fit_mag  fit_std\n');
  for i=1:n_wells,
       jaux=find(w_Id_sh(:)==i);
       
       if (~isempty(jaux))
	 mea_mag = stress_factor*sh_meas_mag(jaux);
	 mea_std = stress_factor*measured.mag_uncertainty(k_sh_mag(jaux));
	 fit_mag = stress_factor*sh_fwd_mag(jaux);
	 fit_std = stress_factor*Ush_fwd_mag(jaux);
	 fprintf('%8.2f %8.2f %8.2f %8.2f\n', [mea_mag(:), mea_std(:), fit_mag(:), fit_std(:)]');
       end
  end
end
   
     
 %%%%%%%%%%%%%%%%%%%% orientation 
 % Sh min azimuth 3 wells raw plots
   handler_shmin_az=figure('Name','Sh azimuth for each well' );
   depth_sh_az=-measured.location(ksh_O,3);   
   w_Id_sh_az=measured.id(ksh_O);
  d_max=max(depth_sh_az);
  d_min=min(depth_sh_az);

  for i=1:n_wells,
      jaux=find(w_Id_sh_az(:)==i);
      
      subplot(1,n_wells,i)
      if (~isempty(jaux))         
       herrorbar(shminaz_mea(jaux),depth_sh_az(jaux),std_shaz_m(jaux),'.r'); hold on;  
       herrorbar(shminaz_hat(jaux),depth_sh_az(jaux),std_shminaz(jaux),'.b'); 
       ylim([d_min d_max])
       set(gca,'yDir','reverse');
       xlim([0 360])
      end 
       title(['Well ID# ', num2str(i)]);
  end
  if ( ~isempty(shminaz_mea) && ~isempty(jaux) )
    %legend(' ','Data',' ', 'Fit ', 'location','Best')
    legend(' ','Measurement Data',' ', 'Inversion Results', 'location','Best')
  else
    fprintf(1, '\n');
    fprintf(2, 'Warning: \n');
    fprintf(1, '	no measurement data for AZ of shmin\n');
    fprintf(1, '	empty plot will be generated\n');
    %fprintf(1, '\n');
  end
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SV plunge....
   handler_sv_pl=figure('Name','Sv Plunge for each well' );
   depth_sv=-measured.location(ksv_O,3);
   w_Id_sv_pl=measured.id(ksv_O);
  d_max=max(depth_sv);
  d_min=min(depth_sv);

  for i=1:n_wells,
      
      jaux=find(w_Id_sv_pl(:)==i);
      
      subplot(1,n_wells,i)
      if (~isempty(jaux))   
       herrorbar(svpl_mea(jaux),depth_sv(jaux),std_svpl_m(jaux),'.r'); hold on;  
       herrorbar(svpl_hat(jaux),depth_sv(jaux),std_svpl(jaux),'.b'); 
       ylim([d_min d_max])
       set(gca,'yDir','reverse');
       xlim([45 135])
      end 
       title(['Well ID# ', num2str(i)]);
  end
  if ( ~isempty(svpl_mea) && ~isempty(jaux) )
    %legend(' ','Data',' ', 'Fit ', 'location','Best')
    legend(' ','Measurement Data',' ', 'Inversion Results', 'location','Best')
  else
    fprintf(1, '\n');
    fprintf(2, 'Warning: \n');
    fprintf(1, '	no measurement data for PL of Sv\n');
    fprintf(1, '	empty plot will be generated\n');
    %fprintf(1, '\n');
  end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
%save plotted figures as fig files
   saveas(handler_sv, 'SvForEachWell.fig', 'fig');
   saveas(handler_shmax, 'ShmaxForEachWell.fig', 'fig');
   saveas(handler_shmin, 'ShminForEachWell.fig', 'fig');
   saveas(handler_shmin_az, 'Shmin_Azimuth_ForEachWell.fig', 'fig');
   saveas(handler_sv_pl, 'Sv_Plunge_ForEachWell.fig', 'fig');

   fprintf(1, '\n');
   fprintf(1, '\n');
   
   return
 

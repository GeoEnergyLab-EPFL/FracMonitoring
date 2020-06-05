function [S_hat]= ComputeUncertainty(matlab_files, mpost, Ctilde, num_ResEleSet,  gauss_number, average_gauss, project_dir)
%
% COMPUTE STRESS UNCERTAINTY
%
% Inputs:
%  matlab_files      :: The filenames of all Matlab format 
%                       stress/strain/displacement saved files
%  mpost             :: 1-dimensional array which contains the inverted 
%                       posterior m
%  Ctilde            :: (m_i, m_j), 2D array which defines the posterior 
%                       covariance matrix
%  num_ResEleSet     :: the number of hex element in the reservoir
%  gauss_number      :: the number of gauss point in each element
%  average_gauss     :: flag if averaging the gauss points or not
%  project_dir       :: directory where the project files are located
%
% Outputs:
%  S_hat             ::(ele_id, stress_id, gauss_point_id), the associated
%                       STD of stress magnitudes
%



   stress_number=6;
   num_model=length(matlab_files);

% Stress_XYZ=reshape(stress_filtered, [num_Stress_read, stress_number, gauss_number]);

%compose sigma_hat by mpost
   num_m=length(mpost);
   if(num_m~=num_model-2)
      error('error: number of mpost is not consistent with Mat results');
   end

   if(average_gauss==1)
      output_gauss_number=1;
   else
      output_gauss_number=gauss_number;
   end

    S_hat=zeros(num_ResEleSet, stress_number, output_gauss_number);
  
    G_aux=zeros(stress_number, num_m);

   for i=1:num_model
      if(file_exist([project_dir, '/', matlab_files{i}.stress_mat])==0)
         error(['One of the matlab_files files does not exist:', matlab_files{i}.stress_mat]);
      end
   end  

   fprintf('Composing G_stress_reservoir\n');
   for i=1:num_model
      fprintf('\t opening Matlab result file: %s \n', matlab_files{i}.stress_mat);
      if(exist('G_stress_reservoir'))
         clear G_stress_reservoir;
      end
      load([project_dir, '/', matlab_files{i}.stress_mat]);

      if(~exist('G_stress_reservoir'))
         error(['G_stress_reservoir not found in ', matlab_files{i}.stress_mat]);
      end
      G_reservoir{i}=G_stress_reservoir;
   end

   %%% loop On Element
   for iele=1:num_ResEleSet
       
       % loop on gauss pt
      for jgauss=1:output_gauss_number
         %compute G_aux for each stress gauss node
         for km=1:num_m
            for lstress=1:stress_number
               G_aux(lstress,km)=G_reservoir{km+2}(iele,lstress,jgauss);
            end
         end
         
           
          %%%% loop on stress component
          for kstress=1:stress_number
              
              cij=G_aux(kstress,:)*Ctilde*G_aux(kstress,:)';
                          
              S_hat(iele,kstress,jgauss)=sqrt(cij);
              
              
          end
       end
   end



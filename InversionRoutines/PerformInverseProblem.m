function [accept,X,PX,k_o,stab]=PerformInverseProblem(mprior,inv_opt)
%
% function where all the steps of the Inverse Problem are processed
% At this stage :: mprior has either been obtained from 1D strain or
% stress model (see EstimateM_strain or EstimateM_stress)
%
%  another possible option could be that it has been directly input from
%  the user as well as the variance (although such a case has still be
%  though of ) so let's leave it for now
%
%
% Define prior covariance matrix (as it s not yet defined)::: a lose
% diagonal matrix scaled by the corresponding mprior values
%
% Inputs:
%  mprior  :: prior m
%  inv_opt :: options used in inversion process
% 
% 
    global Cp ;
    global prior_information measured ;
    global inv_Cd_mag inv_Cd_orient;
    global inv_Cp_sv inv_Cp_sh inv_Cp_shmin_sv_ratio inv_Cp_shmax_sv_ratio;
    global inv_Cp_Q;
    
    fid_log=inv_opt.fid_log;
    
    q=length(mprior);

    stprior=abs(mprior)*100; % "flat" like prior by default
    %stprior = max(abs(mprior))*1000 * ones(size(mprior)); % "flat" like prior by default
  
    %%% check for the case of zero_bending    
    
    % NOTE    shall we do the same before hand to modify mprior after
    % EstimateM_strain and EstimateM_stress ?
    if inv_opt.zero_bending   % over-ride default
        stprior(4:6)=1e-3;  
        mprior(4:6)=0;
    end
    
     Cp=diag(stprior.^2);


% init  inv_Cd_mag,inv_Cd_orient,inv_Cp_sv,inv_Cp_sh,inv_Cp_shmin_sv_ratio,inv_Cp_shmax_sv_ratio (set as globals)
%[inv_Cd_mag,inv_Cd_orient,inv_Cp_sv,inv_Cp_sh,inv_Cp_shmin_sv_ratio,inv_Cp_shmax_sv_ratio]=InitializeCova(measured, prior_information, prior_information.Ns);
[inv_Cd_mag,inv_Cd_orient,inv_Cp_sv,inv_Cp_sh,inv_Cp_shmin_sv_ratio,inv_Cp_shmax_sv_ratio,inv_Cp_Q]=InitializeCova(measured, prior_information, prior_information.Ns);

%%% STEP 1 : Perform a fminsearch from the prior value
fprintf('\nFirst Step of the inversion: a simplex starting from the prior value (1D stress or strain estimate) \n');
fprintf(' --- ' );
fprintf(fid_log, '\nFirst Step of the inversion: a simplex starting from the prior value (1D stress or strain estimate) \n');
fprintf(fid_log, ' --- ' );
%%% SOLUTION VIA A MINIMIZATION
%%%% finding the minimum of Posterior (-log of Posterior) using a
%%%% Nelder-Mead Simplex method :: fminsearchq*200

options = optimset('Display','iter','MaxFunEvals',q*5000,'MaxIter',q*200,'TolFun',1e-3,'TolX',1e-2);
[mpost,fval,exitflag,output]=fminsearch('MinusLogPosterior',mprior,options);


fprintf( '%s\n', output.algorithm);
fprintf('\tfuncCount =  \t\t%d\n', output.funcCount);
fprintf( '\tIterations =  \t\t%d\n', output.iterations);
%fprintf( 'mpost=\n\t\t');
%fprintf('%f\t\t',mpost');
fprintf('\n');
fprintf( 'mpost=\n');
fprintf('\t\t %f\n', mpost');
fprintf('\n');
%fprintf( '\tcurrent mLogP = \n\t\t%f\n', fval);
fprintf( 'current mLogP = \n\t\t%f\n', fval);
fprintf( '%s\n', output.message);

fprintf(fid_log, '%s\n', output.algorithm);
fprintf(fid_log, '\tfuncCount =  \t\t%d\n', output.funcCount);
fprintf(fid_log, '\tIterations =  \t\t%d\n', output.iterations);
%fprintf(fid_log, 'mpost=\n\t\t');
%fprintf(fid_log, '%f\t\t',mpost');
fprintf(fid_log, '\n');
fprintf(fid_log, 'mpost=\n');
fprintf(fid_log, '\t\t %f\n',mpost');
fprintf(fid_log, '\n');
%fprintf(fid_log, '\tcurrent mLogP = \n\t\t%f\n', fval);
fprintf(fid_log, 'current mLogP = \n\t\t%f\n', fval);
fprintf(fid_log, '%s\n', output.message);

%%% STEP 2 : perform a MCMC with automatic stopping criteria in order to
%%% properly sample the posterior pdf

if(inv_opt.MCMC==1)
    % MCMC starting from mpost but with Cp jumps
    mstart=mpost;
    Cstart=Cp;
    plot_flag=0;                % no Plot
    
    % update UserManual accordingly
    %length_default = 15000 ;
    length_default = 10000 ;
    
    fprintf('\n------------------------------------------------ \n' );
    fprintf('2nd Step: MCMC for a robust posterior pdf sampling ' );
    fprintf('\n Markov-Chain-Monte-Carlo on -log( posterior pdf) \n' );
    fprintf(fid_log,'\n------------------------------------------------ \n ' );
    fprintf(fid_log, '2nd Step: MCMC for a robust posterior pdf sampling ' );
    fprintf(fid_log, '\n Markov-Chain-Monte-Carlo on -log( posterior pdf) \n' );
    
    if (inv_opt.MCMC_length==0) %default Automatic arrest of the chain and fixed max length
        Smax=length_default*q; AutomaticStop=1;
        %fprintf('We will perform a MCMC with a maximum number of steps of %i (6 (problem dimension) times 10000) \n',Smax);
        %fprintf(fid_log,'We will perform a MCMC with a maximum number of steps of %i (6 (problem dimension) times 10000)\n',Smax);
        fprintf(1,       'We will perform a MCMC with a maximum number of steps of %i (problem dimension %d times default MCMC_length %d)\n',Smax, q, length_default);
        fprintf(fid_log, 'We will perform a MCMC with a maximum number of steps of %i (problem dimension %d times default MCMC_length %d)\n',Smax, q, length_default);
        
        fprintf(' Heuristic Automatic stopping criteria of the chain is on \n');
        fprintf(fid_log,'Heuristic Automatic stopping criteria of the chain is on  \n');
        
    else
         % Given fixed max length  
        if ( inv_opt.MCMC_length < length_default )
	    length_User = inv_opt.MCMC_length;

            %fprintf(' Given MCMC_length is lower than default value of 10000, we switch back to default.\n');
            %fprintf(fid_log,' Given MCMC_length is lower than default value of 10000, we switch back to default.\n');

            %inv_opt.MCMC_length = length_default;
            %fprintf(        'User-defined MCMC_length (%d) is lower than the default value (%d), we switch back to default.\n', length_User, length_default);
            %fprintf(fid_log,'User-defined MCMC_length (%d) is lower than the default value (%d), we switch back to default.\n', length_User, length_default);

	    fprintf(2, 'Warning: \n');
            fprintf(        'User-defined MCMC_length (%d) is lower than the default value (%d), which might be too small to get a stalbe chain.\n', length_User, length_default);
            fprintf(fid_log,'User-defined MCMC_length (%d) is lower than the default value (%d), which might be too small to get a stalbe chain.\n', length_User, length_default);
        end
        
        Smax=(inv_opt.MCMC_length)*q;AutomaticStop=0; % no arrest control
        fprintf(1,       'We will perform a MCMC with a maximum number of steps of %i (problem dimension %d times user-defined MCMC_length %d) \n',Smax, q, inv_opt.MCMC_length);
        fprintf(fid_log, 'We will perform a MCMC with a maximum number of steps of %i (problem dimension %d times user-defined MCMC_length %d) \n',Smax, q, inv_opt.MCMC_length);
        
        fprintf(' Heuristic Automatic stopping criteria of the chain is disabled \n');
        fprintf(fid_log,'Heuristic Automatic stopping criteria of the chain is disabled  \n');
                
    end


    global info_CP_save;

    if ( isfield(info_CP_save,'MCMC_continue') && info_CP_save.MCMC_continue )
      mcmc_resfile_fullpath = [info_CP_save.project_dir, '/', info_CP_save.mcmc_resfile];
      if ( exist(mcmc_resfile_fullpath, 'file') )
	fprintf(1, '\n');
	fprintf(2, 'Warning: \n');

	question = 'loading previous MCMC results for continuation?';
	auto_skip_answer = 'check';
	IsEnableEmpty = 0;
	flag_YesNo = Smart_AutoAnswer(question, auto_skip_answer, IsEnableEmpty);

	if ( flag_YesNo )
	  load(mcmc_resfile_fullpath, '-MAT');

	  tmp_ndx = find(accept);
	  if ( isempty(tmp_ndx) )
	    fprintf(1, '\n');
	    fprintf(2, 'Warning: \n');
	    fprintf(1, '	no stable result found in %s\n', mcmc_resfile_fullpath);
	    fprintf(1, '	skip loading old MCMC\n');
	  else
	    mstart = X(end,:);
	    if ( isequal(size(mstart), size(mprior)) )
	      ;
	    elseif ( isequal(size(mstart'), size(mprior)) )
	      mstart = mstart';
	    else
	      usr_warn_value(mprior);
	      usr_warn_value(mstart);
	      error('size not match');
	    end
	    fprintf(1, '\n');
	    fprintf(1, 'backup existing MCMC results\n');
	    copyfile(mcmc_resfile_fullpath, [mcmc_resfile_fullpath, '.bak'])
	  end
	end
      end
    end
    fprintf(1, '\n');

    
    tic;
    [accept,X,PX,cf,k_o,stab]=MarkovChainMonteCarlo(fid_log, mstart,Cstart,Smax,plot_flag,AutomaticStop,'MinusLogPosterior');
    toc;
    
    fprintf('\t End of Markov-Chain Monte-Carlo  \n' );
    fprintf('\n------------------------------------------------ \n' );
    fprintf(fid_log, '\t End of Markov-Chain Monte-Carlo  \n' );
    fprintf(fid_log,'\n------------------------------------------------ \n ' );
    if length(PX)<=1
        X=mpost; PX=fval;
    end
    
    
else
    
    fprintf('\t You have chosen NOT to perform a Markov-Chain-Monte-Carlo  \n' );
    fprintf('Be careful, the results (obtained via Simplex algorithm) are most probably not the most probable estimate  \n' );
    fprintf('Uncertainties will be estimated by a simple approximation of the posterior pdf by a Gaussian \n' );
    fprintf(' around the mode obtained by the Simplex algorithm\n');
    fprintf(' for better results (accuracy & robustness) : we recommend to perform a MCMC\n' );
    
    fprintf(fid_log, '\t You have chosen NOT to perform a Markov-Chain-Monte-Carlo  \n' );
    fprintf(fid_log, 'Be careful, the results (obtained via Simplex algorithm) are most probably not the most probable estimate  \n' );
    fprintf(fid_log, 'Uncertainties will be estimated by a simple approximation of the posterior pdf by a Gaussian \n' );
    fprintf(fid_log, ' around the mode obtained by the Simplex algorithm\n');
    fprintf(fid_log, ' for better results (accuracy & robustness) : we recommend to perform a MCMC\n' );
    
    accept=1;
    X=mpost;
    PX=fval;
    k_o=0;
    stab=0;
    
end,


return

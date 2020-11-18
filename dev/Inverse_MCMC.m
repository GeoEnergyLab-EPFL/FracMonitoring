% forward test script 
% Dong Liu -- 29/05/2020

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location
datastor = 'local'; % 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'


%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        %datapath = ['/home/' username(1:end-1) '/data/'];
        %datapath = ['/Users/bricelecampion/Documents/Work/Geomechanics/HydraulicFracturing/Experiments/local_data/'];
        datapath = ['/Users/dongliu/AcousticHF/'];

end

% 2019 acquisitions
datayear = 19;
% test on gabbro
datamonth = 03;
dataday = 14;
starttime = '093621';

% data folder name from experiment date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];

% extract header info from JSON file
fjson = [datapath datafold num2str(starttime) '.json'];
[jsonhdr,myTransducers,myPlattens,myBlock] = load_header(fjson);

%% Read the diffracted arrival
% diffraction data
% we read only the PP diffraction, one can change the wave_type by changing its value
wave_type = 'PP';
sidemarker = ['NV';'SV';'EV';'WV';'NN';'SS';'WW';'EE';'NE';'NW';'SE';'SW';'TT'];
fpath = [datapath datafold ];%'diffraction_picks/'];

% Set basic parameters for inverse problem
gabbro = IsotropicSolid(3000,103.809*1e9,0.284314); % since 28/11/2019 

% set the input parameters
global d  SRPairs ray_type

sig_d = (0.5*1e-6);    % variance of measurement
% this should be the variance of the picked arrival time, can change from case to case
% we adopt here the average picked error 0.5\mu s

global Solid;
Solid = gabbro;

global prior

global m_ind % model indicator

%% M1_Ellipse model indicator as 1
%   a b    x y z, theta, phi psi
% guessed values vector a,b, center coordinate (XYZ),
% (alpha-local rotation; beta-tilt angle, most constrained; gamma-direction of the slope)
m_ind=1;
mp= [ log(0.05);log(0.05); .125;.125;.1285; 0.;0.;0. ];
sig_p = [-log(0.125);-log(0.125);0.02;0.02;0.006;pi/4;pi/60;pi/2]; % guessed variances
% relaxed tilt angle
%sig_p = [.125;.125;0.02;0.02;0.006;pi/4;pi/18;pi/2]; % guessed variances

%% M2_Radial model indicator as 2
% for the radial case r x y z alpha beta
m_ind=2;
mp= [ log(0.05);.125;0.125;.1285;0.;0.];
sig_p = [-log(0.125);0.02;0.02;0.006;pi/60;pi/2]; % guessed variances

%% M3_Ellipse with zero dip model indicator as 3
% for the elliptical case a b x y z self-rotation alpha
m_ind=3;
mp= [ log(0.05);log(0.05); .125;.125;.1285; 0.];
sig_p = [-log(0.125);-log(0.125);0.02;0.02;0.006;pi/4]; % guessed variances

%% M4_Radial with zero dip, model indicator as 4
% for the radial case with fixed z_c coordiante
m_ind=4;
mp= [ log(.05);.125;0.125;.1285];
sig_p = [-log(.125);0.02;0.02;0.006]; % guessed variances

%% M3bis_Ellipse with fixed z-coordinate for center model indicator as 5
% for the radial case r x y z alpha beta
m_ind=5;
global z_c
z_c=0.1285; % measured from the bottom
mp= [ 0.05;0.05;.125;.125;0.; 0.;0.];
sig_p = [.125;0.125;0.02;0.02;pi/4;pi/60;pi/2]; % guessed variances

%% M4bis_Radial with fixed z-coordinate for center model indicator as 6
% for the radial case with fixed z_c coordiante
m_ind=6;
global z_c
z_c=0.1285; % measured from the bottom
mp= [ .05;.125;.125;0.;0.];
sig_p = [.125;0.02;0.02;pi/60;pi/2]; % guessed variances

%% Set the prior
prior = GaussianPrior(mp,sig_p);% build the object

%% Set the sequence, wave_type, block side and its path%
seqnb =  50; % this should be the global sequence number
% Read the SRmap and arrival time
[Pair_info, Pair_acqT]=load_diffraction(fpath, sidemarker, wave_type, seqnb);
SRdiff = SourceReceiverPairs(myTransducers,myPlattens,Pair_info(:,1:2));
d = Pair_info(1:end,3)*1e-6; % arrival time data from diffraction should in s
SRPairs = SRdiff; % SR pairs selected
Cdinvdiag =(1/(1*sig_d^2))*ones(length(d),1); % data inverse of variance
ray_type = ones(size(Pair_info,1),1);

%%  direct unconstrained optimization (Nelder-Mead simplex)
options = optimset('MaxFunEvals',1e5);
[z_opt,bestval,exitflag,output] = fminsearch(@minusPosteriorPDF,[mp;log(sig_d)],options);
% quadratic approximation of the posterior variaance
[Ctilde]=Posterior_Covariance_Matrix_ellipse(z_opt);
sig_app_x=diag(Ctilde).^0.5;
z_opt'
sig_app_x' 

%% DE (differential evolution) algorithm (global optimization)
VTR = -1e6; % stop limit for the objective function, the value of the cost function below which the optimization would be stopped
D = length(mp)+1;% number of parameters of the objective function
           
switch m_ind
    case 1
        XVmax = [log(.250);log(.250);.250;.250;.250;pi/4;pi/2;pi/2; log(20*sig_d)]'; % vector of upper bounds
        XVmin = [log(0.01);log(0.01);0.;0.;0.;-pi/4;-pi/2;-pi/2;log(0.01*sig_d)]';
    case 2
        XVmax = [log(.250);.250;.250;.250;pi/2;pi/2; log(20*sig_d)]';
        XVmin = [log(0.01);0.;0.;0.;-pi/2;-pi/2;log(0.01*sig_d)]';
    case 3
        XVmax = [log(.250);log(.250);.250;.250;.250;pi/2; log(20*sig_d)]'; % vector of upper bounds
        XVmin = [log(0.01);log(0.01);0.;0.;0.;-pi/2;log(0.01*sig_d)]';
    case 4
        XVmax = [log(.250);.250;.250;.250; log(20*sig_d)]';
        XVmin = [log(0.01);0.;0.;0.;log(0.01*sig_d)]';
    case 5
        % needs to write
        %XVmax = [log(.250);log(.250);.250;.250;pi;pi;pi/2; 10*sig_d]'; % vector of upper bounds
        %XVmin = [-5.;-5.;0.;0.;0.;0.;0.;-100]';
    case 6
        % needs to write
        %XVmax = [log(.250);.250;.250;pi;pi; 10*sig_d]';
        %XVmin = [-5.;0.;0.;0.;0.;-100]';
end

NP = 10*D; % number of population members, 10D is by defaut
itermax = 2000; % maximum number of iterations (generations)
F = 0.8; % DE-stepsize F from interval [0, 2]
CR = 0.5; % crossover probability constant from interval [0, 1]
strategy = 7;
refresh = 10; 
[z_opt,bestval,nfeval] = devec3(@minusPosteriorPDF,VTR,D,XVmin,XVmax,[],NP,itermax,F,CR,strategy,refresh);

%% Calculate the covariance matrix

% quadratic approximation of the posterior variance
[Cpost]=Posterior_Covariance_Matrix_ellipse(z_opt);
sig_app_mpost=diag(Cpost).^0.5; % not sure about this.
    
% calculate the Bayes factor
cp=det(diag(sig_p.^2));
% log version
prob_model=-bestval-0.5*log(cp)+0.5*log(det(Cpost));
% nonlog version
%prob_model(i)=exp(-fval)/(cp.^0.5).*((det(Cpost)).^0.5);%/((2*pi).^((length(d)-1)/2));
   
%%%%% check fit
m = z_opt(1:length(z_opt)-1);
noise = z_opt(length(z_opt));
%SRPairs_multiseq{i_seq} = SRPairs;
    switch m_ind
        case 1
            ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 2
            ell = Radial(m(1),m(2:4),m(5),m(6));
        case 3
            ell = Ellipse(m(1),m(2),m(3:5),m(6),0,0);
        case 4
            ell = Radial(m(1),m(2:4),0,0);
        case 5
            ell = Ellipse(m(1),m(2),[m(3:4);z_c],m(5),m(6),m(7));
        case 6
            ell = Radial(m(1),[m(2:3);z_c],m(4),m(5));
        otherwise
            disp('Please check your input vector');
    end
    
res = diffractionForward(Solid,SRPairs,ell,ray_type);% give one the shortest time needed for diffraction

     
% calculate the likelihood
mLikeli=minusLikelihoodEllipseDiff(z_opt);
% better to use the log(likelihood) version
likeli_model=-mLikeli-0.5*length(d)*log(exp(noise)^2)-length(d)/2*log(2*pi);
 
figure
errorbar([1:length(d)]',d*1e6,ones(length(d),1)*exp(noise)*1e6,'b')
hold on
plot([1:length(d)]',res(:,1)*1e6,'*-r'); % the optimized arrival time
xlabel('Source-Receiver Pair Number') % here the label is not clear it is the diffracted SR pick
ylabel('Arrival Time (\mu s)')
legend('real arrival time','calculated arrival time')
hold on
title(['Seq ' num2str(seqnb)])
    
% calculate the correlation matrix
rho_matrix=ones(length(Cpost));
for i_rho=1:length(Cpost)
   for j_rho=1:length(Cpost)
      rho_matrix(i_rho,j_rho)=sqrt(Cpost(i_rho,j_rho).^2/Cpost(i_rho,i_rho)/Cpost(j_rho,j_rho));
   end
end
mpostde=z_opt;
sigpostde=sig_app_mpost;

%% plot the fracture ellipse with diffracted waves
fig_handle = figure('units','normalized','outerposition',[0 0 1 1]);
fig_handle = fractureShape_plot(m,Solid,SRPairs,ray_type,myBlock,myTransducers,myPlattens,fig_handle)
hold on
title(['\fontsize{40}Seq ' num2str(seqnb)])

%% RUNNING A  Marckov-Chain-Monte-Carlo / Metropolis-Hasting algorithm

mstart=z_opt;
Cstart=diag([1./prior.invCpdiag;1/(-log(sig_d))]);% this needs to be checked, this value should be always positive
Smax=50000*length(mstart);
fid_log = fopen('mcmc.txt','w');

[accept, Xf, PXf, cf ,k_o, stab]=MarkovChainMonteCarlo(fid_log, mstart,Cstart,Smax,1,0,@minusPosteriorPDF);

%% save the mcmc results into a json file
MCMCRecord={};
MCMCRecord.mpostde=mpostde;
MCMCRecord.sigpostde=sigpostde;
MCMCRecord.accept=accept;
MCMCRecord.Xf=Xf;
MCMCRecord.PXf=PXf; 
MCMCRecord.cf=cf; % needs to change here
MCMCRecord.ko=k_o; % needs to change here
MCMCRecord.stab=stab; % needs to change here

txtoSave=jsonencode(MCMCRecord);
fname='MCMCRecord.json'; % needs to change here
fid=fopen(fname,'w');
fwrite(fid,txtoSave,'char');
fclose(fid);

%% Read MCMC results from json file
DMCMC = jsondecode(fileread('MCMCRecord.json'));
mpostde=DMCMC.mpostde;
sigpostde=DMCMC.sigpostde;
accept=DMCMC.accept;
Xf=DMCMC.Xf;
PXf=DMCMC.PXf; 
cf=DMCMC.cf; % needs to change here
k_o=DMCMC.ko; % needs to change here
stab=DMCMC.stab;

%% post process - fast way assuming a single gaussian
acceptrate=mean(accept) % should be between 0.2-0.4
%%%%%%%%%%%we add one line here to add the constant terms in
%%%%%%%%%%%minuslogposterior objective function
PXfnew=PXf+0.5*log(det(diag(1./prior.invCpdiag)))+(length(d)+length(prior.invCpdiag))/2*log(2*pi);
n_resampling=2; % sub_sampling of the accepted chain
[mpost,sigpost]=MCMC_hist(accept,PXfnew,Xf,k_o,stab,n_resampling,mpostde,sigpostde);

%% post process - mcmc
% plot the correlation matrix with mcmc fitting results and MCMC
% approximation results.

n_res=2; % subsampling
[MPost,Ctilde_all,W]=ProcessingMCMC(accept,PXf,Xf,k_o,stab,n_res,@minusPosteriorPDF);

%% plot the correlation matrix with mcmc fitting results and DE results

n_res=2; % subsampling
[MPost,Ctilde_all,W]=ProcessingMCMCCompDE(accept,PXf,Xf,k_o,stab,n_res,mpostde,sigpostde,@minusPosteriorPDF);
%% Visualization of the fracture front

m=mpost;
msig=mpost+sigpost;
fig_handle = figure('units','normalized','outerposition',[0 0 1 1]);
fig_handle = fractureShape_plot(m,Solid,SRPairs,ray_type,myBlock,myTransducers,myPlattens,fig_handle)
hold on
fig_handle = fractureShape_plot(msig,Solid,SRPairs,ray_type,myBlock,myTransducers,myPlattens,fig_handle,'r.-')
hold on
title(['\fontsize{40}Seq ' num2str(seqnb)])

%% check fit
m = mpost;
%SRPairs_multiseq{i_seq} = SRPairs;
    switch m_ind
        case 1
            ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 2
            ell = Radial(m(1),m(2:4),m(5),m(6));
        case 3
            ell = Ellipse(m(1),m(2),m(3:5),m(6),0,0);
        case 4
            ell = Radial(m(1),m(2:4),0,0);
        case 5
            ell = Ellipse(m(1),m(2),[m(3:4);z_c],m(5),m(6),m(7));
        case 6
            ell = Radial(m(1),[m(2:3);z_c],m(4),m(5));
        otherwise
            disp('Please check your input vector');
    end
    
res = diffractionForward(Solid,SRPairs,ell,ray_type);% give one the shortest time needed for diffraction

 
figure
errorbar([1:length(d)]',d*1e6,ones(length(d),1)*exp(noise)*1e6,'b')
hold on
plot([1:length(d)]',res(:,1)*1e6,'*-r'); % the optimized arrival time
xlabel('Source-Receiver Pair Number') % here the label is not clear it is the diffracted SR pick
ylabel('Arrival Time (\mu s)')
legend('real arrival time','calculated arrival time')
hold on
title(['Seq ' num2str(seqnb)])










 
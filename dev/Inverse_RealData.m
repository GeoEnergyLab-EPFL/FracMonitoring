% forward test script 
% Dong Liu -- 08/10/2019

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
% Set the sequence, wave_type, block side and its path%
%i_seq = 4; % change this from 1 to seq_n
seqnb =  85;
% we read only the PP diffraction, one can change the wave_type by changing its value
wave_type = 'PP';
sidemarker = ['NV';'SV';'EV';'WV';'SS';'SE';'SW';'NN';'NE';'NW';'EE';'WW';'TT'];
fpath = [datapath datafold ];%'diffraction_picks/'];

% Read the SRmap and arrival time
[Pair_info, AcSeqT]=load_diffraction(fpath, sidemarker, wave_type, seqnb);

% Check the S-R pairs
% build the SRPair object
SRdiff = SourceReceiverPairs(myTransducers,myPlattens,Pair_info(:,1:2));
% plot all the direct traces
fig_b = plotblockwithplattens(myBlock,myPlattens)
fig_handle = plotdirectrays(SRdiff,fig_b);

%% Set basic parameters for inverse problem
% gabbro = IsotropicSolid(3050,97.5867*1e9,0.3119); % vel in m/s before 28/11/2019
% relation between E, density and Poisson's ratio
gabbro = IsotropicSolid(3000,103.809*1e9,0.284314); % since 28/11/2019 


% set the input parameters
global d  SRPairs ray_type
d = Pair_info(1:end,3)*1e-6; % arrival time data from diffraction should in s
SRPairs = SRdiff; % SR pairs selected
ray_type = ones(size(Pair_info,1),1);

sig_d = (0.5*1e-6);    % variance of measurement
% this should be the variance of the picked arrival time, can change from case to case
% we adopt here the average picked error 0.5\mu s

global Solid  Cdinvdiag;
Solid = gabbro;
Cdinvdiag =(1/(1*sig_d^2))*ones(length(d),1); % data inverse of variance

global prior

global m_ind;

%% M1_Ellipse model indicator as 1
%   a b    x y z, theta, phi psi
% guessed values vector a,b, center coordinate (XYZ),
% (alpha-local rotation; beta-tilt angle, most constrained; gamma-direction of the slope)
m_ind=1;
mp= [ .05;.05; .125;.125;.1285; 0.;0.;0. ];
sig_p = [.125;.125;0.02;0.02;0.006;pi/4;pi/60;pi/2]; % guessed variances

%% M2_Radial model indicator as 2
% for the radial case r x y z alpha beta
m_ind=2;
mp= [ .05;.125;0.125;.1285;0.;0.];
sig_p = [.125;0.02;0.02;0.006;pi/60;pi/2]; % guessed variances

%% M3_Ellipse model indicator as 3
%   a b    x y z, theta, phi psi
% guessed values vector a,b, center coordinate (XYZ),
% (alpha-local rotation; beta-tilt angle, most constrained; gamma-direction of the slope)
m_ind=3;
mp= [ .05;.05; .125;.125;.1285; 0.];
sig_p = [.125;.125;0.02;0.02;0.006;pi/4]; % guessed variances

%% M2_Radial model indicator as 4
% for the radial case r x y z alpha beta
m_ind=4;
mp= [ .05;.125;0.125;.1285];
sig_p = [.125;0.02;0.02;0.006]; % guessed variances

%% M3bis_Ellipse with fixed z-coordinate for center
% for the radial case r x y z alpha beta
m_ind=5;
global z_c
z_c=0.1285; % measured from the bottom
mp= [ 0.05;0.05;.125;.125;0.; 0.;0.];
sig_p = [.125;0.125;0.02;0.02;pi/4;pi/60;pi/2]; % guessed variances

%% M4bis_Radial with fixed z-coordinate for center
% for the radial case with fixed z_c coordiante
m_ind=6;
global z_c
z_c=0.1285; % measured from the bottom
mp= [ .05;.125;.125;0.;0.];
sig_p = [.125;0.02;0.02;pi/60;pi/2]; % guessed variances

%% Set the prior
prior = GaussianPrior(mp,sig_p);% build the object

%%  direct unconstrained optimization (Nelder-Mead simplex)

[x,fval,exitflag,output] = fminsearch(@minusPosteriorPDF,mp);
% quadratic approximation of the posterior variaance
[Ctilde]=Posterior_Covariance_Matrix_ellipse(x);
sig_app_x=diag(Ctilde).^0.5;
x'
sig_app_x' 

%% DE (differential evolution) algorithm (global optimization)
VTR = 0.1; % stop limit for the objective function, the value of the cost function below which the optimization would be stopped
D = length(mp);% number of parameters of the objective function
           
XVmin = 0.*mp'; % vector of lower bounds XVmin(1) ... XVmin(D)
switch m_ind
    case 1
        XVmax = [.250;.250;.250;.250;.250;pi;pi;pi/2]'; % vector of upper bounds
    case 2
        XVmax = [.250;.250;.250;.250;pi;pi]';
    case 3
        XVmax = [.250;.250;.250;.250;.250;pi]'; % vector of upper bounds
    case 4
        XVmax = [.250;.250;.250;.250]';    
    case 5
        XVmax = [.250;.250;.250;.250;pi;pi;pi/2]'; % vector of upper bounds
    case 6
        XVmax = [.250;.250;.250;pi;pi]';
end

NP = 10*D; % number of population members, 10D is by defaut
itermax = 2000; % maximum number of iterations (generations)
F = 0.8; % DE-stepsize F from interval [0, 2]
CR = 0.5; % crossover probability constant from interval [0, 1]
strategy = 7;
refresh = 10; 
[bestmem,bestval,nfeval] = devec3(@minusPosteriorPDF,VTR,D,XVmin,XVmax,[],NP,itermax,F,CR,strategy,refresh);

%% compare the arrival time from the optmization and the measurement
m = bestmem';
switch m_ind
    case 1 % ellipse case
        ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
    case 2 % radial case
        ell = Radial(m(1),m(2:4),m(5),m(6));
    case 3 % ellipse case
        ell = Ellipse(m(1),m(2),m(3:5),m(6),0,0);
    case 4 % radial case
        ell = Radial(m(1),m(2:4),0,0);
    case 5
        ell = Ellipse(m(1),m(2),[m(3:4),z_c],m(5),m(6),m(7));
    case 6
        ell = Radial(m(1),[m(2:3),z_c],m(4),m(5));
end 
res = diffractionForward(Solid,SRPairs,ell,ray_type);% give one the shortest time needed for diffraction
figure
hold on
title('model vs data');
%plot(res(:,1)*1e6,'ob')
errorbar(res(:,1)*1e6,d*1e6,ones(length(d),1)*sig_d*1e6,'b')

figure
errorbar([1:length(d)]',d*1e6,ones(length(d),1)*sig_d*1e6,'b')
hold on
plot([1:length(d)]',res(:,1)*1e6,'*-r'); % the optimized arrival time
xlabel('Different Source-Receiver Pairs') 
ylabel('Arrival Time (\mu s)')
legend('real arrival time','calculated arrival time')


%% plot the fracture ellipse with diffracted waves
fig_handle = figure('units','normalized','outerposition',[0 0 1 1]);
fig_handle = fractureShape_plot(m,Solid,SRPairs,ray_type,myBlock,myTransducers,myPlattens,fig_handle)
hold on
title(['\fontsize{40}Seq ' num2str(seqnb) ': ' datestr(AcSeqT)])

%% RUNNING A  Marckov-Chain-Monte-Carlo / Metropolis-Hasting algorithm

mstart=m;
Cstart=diag(1./prior.invCpdiag);
Smax=50000*length(m);
fid_log = fopen('mcmc.txt','w');

[accept, Xf, PXf, cf ,k_o, stab]=MarkovChainMonteCarlo(fid_log, mstart,Cstart,Smax,1,0,@minusPosteriorPDF);


%% post process - fast way assuming a single gaussian
acceptrate=mean(accept) % should be between 0.2-0.4

n_resampling=1; % sub_sampling of the accepted chain
[mpost,sigpost]=MCMC_hist(accept,PXf,Xf,k_o,stab,n_resampling);

%% post process - mcmc

%strings_m={'a','b','x_c','y_c','z_c','phi','theta','psi'};
%col={'.r','.b','.m','.g','.y','.k','+r','+b'};

%[MPost,Ctilde_all,W]=fitGaussiansMixture(Xsp,@minusPosteriorPDF);
n_res=2; % subsampling
[MPost,Ctilde_all,W]=ProcessingMCMC(accept,PXf,Xf,k_o,stab,n_res,@minusPosteriorPDF);

%% Comparison betweeen the models
Nd=length(d);
k_s=[k_o:n_res:length(PXf)];
Xsp=Xf(k_s,:);
n_sample_ellipse=length(Xsp(:,1));
P_ellipse=(1/(1*sig_d))*exp(-minusLikelihoodEllipseDiff(MPost))*(1/2/pi)^(Nd/2);
%%
n_res=2; % subsampling
k_s=[k_o:n_res:length(PXf)];
Xsp=Xf(k_s,:);
n_sample_radial=length(Xsp(:,1));
[MPost_radial,Ctilde_all,W]=ProcessingMCMC(accept,PXf,Xf,k_o,stab,n_res,@minusPosteriorPDF);
P_radial=(1/(1*sig_d))*exp(-minusLikelihoodEllipseDiff(MPost_radial))*(1/2/pi)^(Nd/2);
BIC_radial=6*log(n_sample_radial)-2*log(P_radial);
BIC_ellipse=8*log(n_sample_ellipse)-2*log(P_ellipse);
%%
% quadratic approximation of the posterior variance
[Ctilde]=Posterior_Covariance_Matrix_ellipse(mpost);
sig_app=diag(Ctilde).^0.5;
sig_app'
sigpost'

% compare the best fit from the differential evolution
[Ctilde]=Posterior_Covariance_Matrix_ellipse(bestmem');
sig_app=diag(Ctilde).^0.5;

sig_app'
sigpost'
sig_app_x'
%%
 
m=mpost
msig=mpost+sigpost;

fig_handle = figure('units','normalized','outerposition',[0 0 1 1]);
fig_handle = fractureShape_plot(m,Solid,SRPairs,ray_type,myBlock,myTransducers,myPlattens,fig_handle)
hold on
fig_handle = fractureShape_plot(msig,Solid,SRPairs,ray_type,myBlock,myTransducers,myPlattens,fig_handle,'r.-')
hold on
title(['\fontsize{40}Seq ' num2str(seqnb) ': ' datestr(AcSeqT)])

%%  saving stuff from dong - to be improved

% Initialize the saved file
res_multiseq = cell(seq_n,5); % sequence number + m-vector 

%% Arrange the results into array or cell and begin the calculation of the next sequence
res_multiseq{i_seq,1} = seqnb; 
res_multiseq{i_seq,2} = AcSeqT;
res_multiseq{i_seq,3} = SRPairs;
res_multiseq{i_seq,4} = m';
res_multiseq{i_seq,5} = m';% m-vector from Monte-Carlo method
%% Save the results into a json file
% One can save the data after running all the sequences one wants to look at
% gloabl sequence  number
for j=1:seq_n
DiffRecord(j).seqnb=res_multiseq{j,1};
% solid properties
DiffRecord(j).solid=gabbro;
% acquisition time
DiffRecord(j).acqT=res_multiseq{j,2};
% SR_map used
DiffRecord(j).SRmap=res_multiseq{j,3}.SRmap;
DiffRecord(j).mDE=res_multiseq{j,4};
DiffRecord(j).mMCMC=res_multiseq{j,5};% res_multiseq{j,5} maybe?
DiffRecord(j).pickedArrival=Pair_info(1:end,3)*1e-6;
% Simulation parameters
DiffRecord(j).variance=vv;
DiffRecord(j).mtrial=mp;
DiffRecord(j).sigmatrial=sig_p;
end
% write into the json file 
txtoSave=jsonencode(DiffRecord);
fname='testJsonSave.json';
fid=fopen(fname,'w');
fwrite(fid,txtoSave,'char');
fclose(fid);

%% Read the results from the saved json file
filename = 'testJsonSave.json';
Shape_fracture = jsondecode(fileread(filename));
nb_seq=length(Shape_fracture);

Material_gabbro = IsotropicSolid(Shape_fracture(1).solid.density,...
    Shape_fracture(1).solid.Young,Shape_fracture(1).solid.nu);

% plot the ellipse fracture and diffracted waves
fig2 = figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:nb_seq
    m_i = Shape_fracture(i).mDE;
    SRdiff_i = SourceReceiverPairs(myTransducers,myPlattens,Shape_fracture(i).SRmap);
    fractureShape_plot(m_i,Material_gabbro,SRdiff_i,ones(size(Shape_fracture(i).SRmap,1),1),...
        myBlock,myTransducers,myPlattens,fig2);%,[],widthseqpath, widthopeningpath);
    hold on
    title(['\fontsize{40}Seq ' num2str(Shape_fracture(i).seqnb) ': ' Shape_fracture(i).acqT])
    pause(2)
    if i<nb_seq
        clf;
    end
end

%% plot with the fracture opening profile from the transmission
widthPath = '/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/G01SCCERWidth/'; % we load the data from the SCCER-Soe Conference
widthseqpath = [widthPath 'OpeningSequenceSCCERG01.txt'];
widthopeningpath = [widthPath 'OpeningSCCERG01.txt'];

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
fig_handle = fractureShape_plot(m,Solid,SRPairs,ray_type,myBlock,myTransducers,myPlattens,fig1,[],widthseqpath,widthopeningpath,seqnb);
hold on
title(['\fontsize{40}Seq ' num2str(seqnb) ': ' datestr(AcSeqT)])



 
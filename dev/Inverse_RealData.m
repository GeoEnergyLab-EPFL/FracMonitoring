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
        datapath = ['/Users/bricelecampion/Documents/Work/Geomechanics/HydraulicFracturing/Experiments/local_data/'];

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

% load timing data
%fbin = [datapath datafold num2str(starttime) '.bin'];
%AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format

% %% Set the sequences that you would like to look at
% seq = [30 50 70 90];
% seq_n = length(seq);
% res_multiseq = zeros(seq_n,9); % sequence number + m-vector 
% SRPairs_multiseq = cell(seq_n,1);
% % time of the loaded sequence
% %AcSeqT = AcqTime(seq);

%% Read the diffracted arrival
% Set the sequence, wave_type, block side and its path%
%i_seq = 4; % change this from 1 to seq_n
seqnb =  60;
% we read only the PP diffraction, one can change the wave_type by changing its value
wave_type = 'PP';
sidemarker = ['N','S','E','W'];
fpath = [datapath datafold 'diffraction_picks/'];

% Read the SRmap and arrival time
Pair_info=load_diffraction(fpath, sidemarker, wave_type, seqnb);

% Check the S-R pairs
% build the SRPair object
SRdiff = SourceReceiverPairs(myTransducers,myPlattens,Pair_info(:,1:2));
% plot all the direct traces
fig_b = plotblockwithplattens(myBlock,myPlattens)
fig_handle = plotdirectrays(SRdiff,fig_b);

%% Set basic parameters for inverse problem
gabbro = IsotropicSolid(3050,97.5867*1e9,0.3119); % vel in m/s
% relation between E, density and Poisson's ratio

% set the input parameters
global d  SRPairs 
d = Pair_info(1:end,3)*1e-6; % arrival time data from diffraction should in s
SRPairs = SRdiff; % SR pairs selected

sig_d = (0.5*1e-6);    % variance of measurement
% this should be the variance of the picked arrival time, can change from case to case
% we adopt here the average picked error 0.5\mu s

global Solid  Cdinvdiag;
Solid = gabbro;
Cdinvdiag =(1/(1*sig_d^2))*ones(length(d),1); % data inverse of variance

global prior

%   a b    x y z, theta, phi psi 
% guessed values vector a,b, center coordinate (XYZ), euler angle (alpha beta gamma)
mp= [ .05;.05; .125;.125;.125; 0.;0.;0. ]; 
sig_p = [.125;.125;0.02;0.02;0.005;pi/4;pi/40;pi/40]; % guessed variances


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
XVmax = [.250;.250;.250;.250;.250;pi/2;pi/2;pi]'; % vector of upper bounds XVmax(1) ... XVmax(D)
NP = 10*D; % number of population members, 10D is by defaut
itermax = 2000; % maximum number of iterations (generations)
F = 0.8; % DE-stepsize F from interval [0, 2]
CR = 0.5; % crossover probability constant from interval [0, 1]
strategy = 7;
refresh = 10; 
[bestmem,bestval,nfeval] = devec3(@minusPosteriorPDF,VTR,D,XVmin,XVmax,[],NP,itermax,F,CR,strategy,refresh);

%% compare the arrival time from the optmization and the measurement
m = bestmem';
%res_multiseq(i_seq,1:end) = [seqnb m'];

%SRPairs_multiseq{i_seq} = SRPairs;
ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
res = diffractionForward(Solid,SRPairs,ell);% give one the shortest time needed for diffraction
figure
title('model vs data');
%plot(res(:,1)*1e6,'ob')
errorbar(res(:,1)*1e6,d*1e6,ones(length(d),1)*sig_d*1e6,'b')

figure
errorbar([1:length(d)]',d*1e6,ones(length(d),1)*sig_d*1e6,'b')
hold on
%plot(d*1e6); hold on; % the real arrival time
%errorbar([1:length(d)]',res(:,1)*1e6,ones(length(d),1)*sig_d*1e6,'b')
plot([1:length(d)]',res(:,1)*1e6,'*-r'); % the optimized arrival time
xlabel('Source-Receiver Pair Number') % here the label is not clear it is the diffracted SR pick
ylabel('Arrival Time (\mu s)')
legend('real arrival time','calculated arrival time')


%%
ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
res = diffractionForward(Solid,SRPairs,ell);% give one the shortest time needed for diffraction

fig_handle = figure;
plotblockwithplattens(myBlock,myPlattens,fig_handle)

hold on
plotdiffrays(SRPairs,res(:,2),res(:,3),res(:,4),fig_handle);% one should plot the trays with the corresponding diffracted points

hold on
plotEllipse(ell,fig_handle,'b.-');
hold on
plot3(res(:,2),res(:,3),res(:,4),'.g','MarkerSize',30);



%% RUNNING A  Marckov-Chain-Monte-Carlo / Metropolis-Hasting algorithm

mstart=m;
Cstart=diag(1./prior.invCpdiag);
Smax=20000*length(m);
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

%%
% quadratic approximation of the posterior variaance
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
ellpost=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
ellpostsig=Ellipse(msig(1),msig(2),msig(3:5),msig(6),msig(7),msig(8));
res= diffractionForward(Solid,SRPairs,ell);

fig_b=plotblockwithplattens(myBlock,myPlattens)

%
fig_handle=plotdirectrays(SRdiff,fig_b);
 
fig_handle=plotEllipse(ellpost,fig_handle,'b');
fig_handle=plotEllipse(ellpostsig,fig_handle,'r-');

hold on
plot3(res(:,2),res(:,3),res(:,4),'.g','MarkerSize',30);


%%  saving stuff from dong - to be improved
%%%%%%%%%%%
%% plot with the fracture opening profile from the transmission
widthPath = '/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/G01SCCERWidth/'; % we load the data from the SCCER-Soe Conference
widthseqpath = [widthPath 'OpeningSequenceSCCERG01.txt'];
widthopeningpath = [widthPath 'OpeningSCCERG01.txt'];
fig1 = figure('units','normalized','outerposition',[0 0 1 1])
% without opening data, one can remove the last two arguments
fig_handle = fractureShape_plot(seqnb,m,Solid,SRPairs,myBlock,myTransducers,myPlattens,AcqTime,fig1, widthseqpath, widthopeningpath)









%% Save the results into a json file
% One can save the data after running all the sequences one wants to look at
% gloabl sequence  number
for j=1:seq_n
DiffRecord(j).seqnb=seq(j);
% solid properties
DiffRecord(j).solid=gabbro;
% acquisition time
DiffRecord(j).acqT=AcSeqT(j);
% SR_map used
DiffRecord(j).SRmap=SRPairs_multiseq{j}.SRmap;
DiffRecord(j).mDE=res_multiseq(j,2:end);
DiffRecord(j).mMCMC=res_multiseq(j,2:end);% resMC
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
    fractureShape_plot(Shape_fracture(i).seqnb,m_i,Material_gabbro,SRdiff_i,...
        myBlock,myTransducers,myPlattens,AcqTime,fig2);%,widthseqpath, widthopeningpath);
    
    pause(2)
    if i<nb_seq
        clf;
    end
end


%%

%%

%plot(Xsp(:,1),Xsp(:,2),'.')

plot(Xsp(:,8),Xsp(:,7),'.')


%%
%%[MPost,Ctilde_all]=ProcessingMCMC(accept,PX,X,k_o,stab,mcmc_res_logfile)

% [MPost,Ctilde_all]=ProcessingMCMC(accept,PXf,Xf,k_o,stab,@minusPosteriorPDF,[]);
 
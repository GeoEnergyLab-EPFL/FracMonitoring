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
seqnb =  79;%85;%40;%79;
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

global Solid
Solid = gabbro;

global prior

global m_ind;

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
sig_p = [-log(0.125);0.02;0.02;0.006;pi/18;pi/2]; % guessed variances

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

[x,fval,exitflag,output] = fminsearch(@minusPosteriorPDF,[mp;sig_d]);
% quadratic approximation of the posterior variaance
[Ctilde]=Posterior_Covariance_Matrix_ellipse(x);
sig_app_x=diag(Ctilde).^0.5;
x'
sig_app_x' 

%% DE (differential evolution) algorithm (global optimization)
VTR = -1e6; % stop limit for the objective function, the value of the cost function below which the optimization would be stopped
D = length(mp)+1;% number of parameters of the objective function
           
%XVmin =[0.*mp' 0]; % vector of lower bounds XVmin(1) ... XVmin(D)
switch m_ind
    case 1
        XVmax = [log(.250);log(.250);.250;.250;.250;pi/2;pi;pi; log(10*sig_d)]'; % vector of upper bounds
        XVmin = [log(0.01);log(0.01);0.;0.;0.;0;0.;0.;log(0.01*sig_d)]';
    case 2
        XVmax = [log(.250);.250;.250;.250;pi;pi; 10*sig_d]';
        XVmin = [log(0.01);0.;0.;0.;0.;0.;-100]';
    case 3
        XVmax = [log(.250);log(.250);.250;.250;.250;pi; 10*sig_d]'; % vector of upper bounds
        XVmin = [-5.;-5.;0.;0.;0.;0.;-100]';
    case 4
        XVmax = [log(.250);.250;.250;.250; 10*sig_d]';
        XVmin = [-5.;0.;0.;0.;-100]';
    case 5
        XVmax = [log(.250);log(.250);.250;.250;pi;pi;pi/2; 10*sig_d]'; % vector of upper bounds
        XVmin = [-5.;-5.;0.;0.;0.;0.;0.;-100]';
    case 6
        XVmax = [log(.250);.250;.250;pi;pi; 10*sig_d]';
        XVmin = [-5.;0.;0.;0.;0.;-100]';
end

NP = 10*D; % number of population members, 10D is by defaut
itermax = 2000; % maximum number of iterations (generations)
F = 0.8; % DE-stepsize F from interval [0, 2]
CR = 0.5; % crossover probability constant from interval [0, 1]
strategy = 7;
refresh = 10; 
[z_opt,bestval,nfeval] = devec3(@minusPosteriorPDF,VTR,D,XVmin,XVmax,[],NP,itermax,F,CR,strategy,refresh);

%% Calculate the covariance matrix

[Cpost]=Posterior_Covariance_Matrix_ellipse(z_opt);
sig_app_mpost=diag(Cpost).^0.5; 

% calculate the Bayes factor
cp=det(diag(sig_p.^2));
% log version
prob_model=-bestval-0.5*cp+0.5*det(Cpost);
% nonlog version
%prob_model(i)=exp(-fval)/(cp.^0.5).*((det(Cpost)).^0.5);%/((2*pi).^((length(d)-1)/2));
   
%%%%% check fit
m = z_opt(1:length(z_opt)-1);
noise = z_opt(length(z_opt));



%% compare the arrival time from the optmization and the measurement
bestmem=z_opt(1:length(z_opt)-1);
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
Smax=40000*length(m);%40000 for Seq40 and Seq79
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

%% Calculate the Person's correlation coefficient
fid=fopen('MCMCmSeq79M1.txt');
m_all=fscanf(fid,'%f',[8 131863]);% the value 986 can be different here
% [8 986] for Seq85
% [8 19150] for Seq40
% [8 131863] for Seq79
%%
corr=corrcoef(m_all');
R=corrplot(m_all');

%%
h=heatmap(abs(corr))
colormap(viridis());

%% plot the colormap for the correlation coefficients
abscorr=abs(corr);

fig_handle = figure('DefaultAxesFontSize',24);
set(gcf,'Position',[100 100 600 600]);
figure(fig_handle)
imagesc(1:size(abscorr,2),1:size(abscorr,1),abscorr)
axis square
caxis([0 1])
colormap(viridis());
colorbar
set(gca,'xtick',[]);
set(gca,'ytick',[]);


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



 
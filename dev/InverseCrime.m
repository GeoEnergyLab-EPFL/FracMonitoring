% forward test script 
% Dong Liu -- 08/10/2019

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location
datastor = 'gel-nas1'; % 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'

%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        datapath = ['/home/' username(1:end-1) '/data/'];
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
fbin = [datapath datafold num2str(starttime) '.bin'];
% AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format

%% Set the sequences that you would like to look at
seq = [30 50 70 90];
seq_n = length(seq);

%% Read the diffracted arrival
% Set the sequence, wave_type, block side and its path
i_seq = 3; % change this from 1 to seq_n
seqnb = seq(i_seq);
% we read only the PP diffraction, one can change the wave_type by changing its value
wave_type = 'PP';
sidemarker = ['N','S','E','W'];
fpath = ['/Users/dongliu/AcousticHF/' datafold];

% Read the SRmap and arrival time
[Pair_info, AcSeqT]=load_diffraction(fpath, sidemarker, wave_type, seqnb);

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

vv = (0.5*1e-6)^2;
% this should be the variance of the picked arrival time, can change from case to case
% we adopt here the average picked error 0.5\mu s

global Solid  Cdinvdiag;
Solid = gabbro;
Cdinvdiag = 1/(1*vv).*(0.*d+1); % data inverse of variance

global prior

mp= [ .05;.05; .125;.125;.125; 0.;0.;0. ]; % guessed values vector a,b, center coordinate (XYZ), euler angle (alpha beta gamma)
sig_p = [.125;.125; 0.02;0.02;0.005;pi/4;pi/40;pi/40]; % guessed variances
prior = GaussianPrior(mp,sig_p);% build the object

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
ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
res = diffractionForward(Solid,SRPairs,ell);% give one the shortest time needed for diffraction

figure
plot(d*1e6); hold on; % the real arrival time
plot(res(:,1)*1e6); % the optimized arrival time
xlabel('Source-Receiver Pair Number')
ylabel('Arrival Time (\mu s)')
legend('real arrival time','calculated arrival time')

%% plot with the fracture opening profile from the transmission
% widthPath = '/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/G01SCCERWidth/'; % we load the data from the SCCER-Soe Conference
% widthseqpath = [widthPath 'OpeningSequenceSCCERG01.txt'];
% widthopeningpath = [widthPath 'OpeningSCCERG01.txt'];
fig1 = figure('units','normalized','outerposition',[0 0 1 1])
% without opening data, one can remove the last three arguments
fig_handle = fractureShape_plot(m,Solid,SRPairs,myBlock,myTransducers,myPlattens, fig1)%, widthseqpath, widthopeningpath, seqnb)
hold on
title(['\fontsize{40}Seq ' num2str(seqnb) ': ' datestr(AcSeqT)])

%% Initialize the saved file
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
    fractureShape_plot(m_i,Material_gabbro,SRdiff_i,...
        myBlock,myTransducers,myPlattens,fig2);%,widthseqpath, widthopeningpath);
    hold on
    title(['\fontsize{40}Seq ' num2str(Shape_fracture(i).seqnb) ': ' Shape_fracture(i).acqT])
    pause(2)
    if i<nb_seq
        clf;
    end
end

%% --------------------------Below is the original script-----------------------------
% Not fully coded yet

%% MCMC -- Markov Chain Monte Carlo Method
% we do not change much this part of the code, please refer to the
% orignal code kept below
mstart=m;
Cstart=diag(1./prior.invCpdiag);
Smax=20000*length(m);
fid_log = fopen('test.txt','w');% record the history

[accept, Xf, PXf, cf ,k_o, stab]=MarkovChainMonteCarlo(fid_log, mstart,Cstart,Smax,1,0,@minusPosteriorPDF);

%% 
acceptrate=mean(accept) % should be between 0.2-0.4;

% resampling step
ns=1;
q=length(Xf(1,:)); % problem dimension


tmp_ndx = find(accept); % non-zero index

k_s=[k_o:ns:length(PXf)];% index where the distribution is stable
  
mLogP=PXf(k_s);
Xsp=Xf(k_s,:);

% -LOG POST histogram
handler_histpost=figure('Name','Histogram of - Log Posterior');
hist(mLogP,60);  % it should look like a lognormal pdf % kind of a probability plot
title('Histogram of - Log Posterior');
legend('It should look like a log Normal pdf' );

%%
figure % plot of the values of a and b, all combinaisions
plot(Xsp(:,1),Xsp(:,2),'.')
xlabel('a (m)')
ylabel('b (m)')

figure % plot of the values of x and y, all combinaisions
plot(Xsp(:,3),Xsp(:,4),'.')
xlabel('x (m)')
ylabel('y (m)')

figure % plot of the Eulerian angles
plot(Xsp(:,8),Xsp(:,7),'.')

%%

%plot(Xsp(:,1),Xsp(:,2),'.')

plot(Xsp(:,8),Xsp(:,7),'.')

%%
%%[MPost,Ctilde_all]=ProcessingMCMC(accept,PX,X,k_o,stab,mcmc_res_logfile)

% [MPost,Ctilde_all]=ProcessingMCMC(accept,PXf,Xf,k_o,stab,@minusPosteriorPDF,[]);
 

%%
inv_opt.MCMC=1;
[mpost,Ctilde] = PostProcessInversion(inv_opt,accept,PXf,Xf,k_o,stab,@minusPosteriorPDF,'test.txt');

%% Plot the results

m_MC=mean(Xsp);

m_MC=mean(Xf(20000:30000,:)); 

ell_MC = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
resMC= diffractionForward(Solid,SRPairs,ell_MC); % calculated arrival time


figure
plot(d*1e6); hold on; % the real arrival time
plot(resMC(:,1)*1e6); % the optimized arrival time
xlabel('Source-Receiver Pair Number')
ylabel('Arrival Time (\mu s)')
legend('real arrival time','calculated arrival time')

% plot the fracture geometry and diffracted waves
fig1 = figure('units','normalized','outerposition',[0 0 1 1])
fig_handle = fractureShape_plot(m,Solid,SRPairs,myBlock,myTransducers,myPlattens,fig1, widthseqpath, widthopeningpath, seqnb)
hold on
title(['\fontsize{40}Seq ' num2str(seqnb) ': ' datestr(AcSeqT)])


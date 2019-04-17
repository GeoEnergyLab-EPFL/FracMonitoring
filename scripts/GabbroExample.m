%% This is an example of a way to load and look at acoustic data from a HF block experiment 
% GEL - EPFL - 2019
% 
%  it will be improved / corrected as features evolve 
%
%
% make sure the src and subfolders are in the matlab path

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location 
% note here use 'gel-nas1' -  
datastor = 'gel-nas1'; % 'local'; 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'

%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        datapath = ['/Users/' username(1:end-1) '/data/']; % here /Users is on OS X, on linux replace by home
end

% 2019 acquisitions

datayear = 19;
%
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


% read active acoustic parameters from JSON header
ActiveAcoustic = jsonhdr.ActiveAcousticInfos;
% extract relevant parameters
np = ActiveAcoustic.NumberOfPoints;
ns = ActiveAcoustic.NumberOfSources;
nr = ActiveAcoustic.NumberOfReceivers;
Fs = ActiveAcoustic.SamplingFrequency_MHz_*1E6;
dt = 1/Fs;  % time step
t0 = 0;     % initial time
T = t0+dt*linspace(0,np-1,np)'; % time vector
Fn = 0.5*Fs;    % Nyquist frequency (Hz)



%%  read acoustic data from tge binary file 
fbin = [datapath datafold num2str(starttime) '.bin'];
% note it is good practice not to load all the data  / i.e. all the sequences
% 1 sequence has np \times 32 \times 32 data points => 1024*np*8 Bytes (64
% bit storage) 
% default param np=8000 so about 6.5Mb for 1 sequence of acquistion !

[dataseq1, nq] = load_data(fbin,1:2:60,1:32,1:32);  % here we load  30 sequences from 1 to 60 by 2

%% load timing data
% this load the time of each acoustic acquisition sequence 
% note it is the time at which sequence ends and is stored to the binary
% file

AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format
% note here this load for all sequence -> API could be change for case like
% this one where we do not load all sequence

% time of the loaded sequence
AcSeqT=AcqTime(1:2:60);

tt=datenum(AcSeqT)-datenum(AcSeqT(1))
%% loadin the low rate data 
[lowrateData, lowrateTime, lowrateHeader] = load_lowrate(fbin);

% plot injection flows
figure
plot(lowrateTime,lowrateData(:,11),lowrateTime,lowrateData(:,14))
axis tight
xlim([lowrateTime(1) lowrateTime(end)])
ylim([0 1])

% plot pressure gauge and active acq times
figure
plot(lowrateTime,lowrateData(:,2),lowrateTime,lowrateData(:,3))
axis tight
xlim([lowrateTime(1) lowrateTime(end)])
ylim([0 40])
% add acq times
hold on
plot(AcqTime,ones(size(AcqTime))*10,'k.')
hold on
plot(AcSeqT,ones(size(AcSeqT))*10,'r.')


% for better plot we need time in seconds and define an initial zero time 

% here plot of flat-jacks pairs ... to help the analysis.


%% Quality Control of acoustic signal: SNR matrix on all pairs.

endnoise = 100; % burning period for EM noise from the source excitation
[normAll, h_SNR]  = signal_strength(squeeze(dataseq1(1,:,:,:)),endnoise)


%%  Selection of source - receiver pairs to compute velocities at t=0

% We choose all pairs from opposite platen 
oppSR=AllPairsOppositePlattens(myTransducers,myPlattens);

fig_b=plotblockwithplattens(myBlock,myPlattens);
%
fig_handle=plotdirectrays(oppSR,fig_b);

%% get the S-R pairs  with sufficient snr

snrv=zeros(length(oppSR.SRmap),1); % here I do a loop because Matlab is stupid / not well designed for array access
for i=1:length(snrv)
    snrv(i)=normAll(oppSR.SRmap(i,1),oppSR.SRmap(i,2));
end

figure 
plot(snrv)

[r c ~]=find(snrv>2);  % adjust threshold of snr above which we look at the data - we are strict here

snrvf=snrv(r);

% create a new object source-receiver pairs containing the pairs from
% opposite platen with a sufficient SNR 
submap= oppSR.SRmap(r,:);
subwave=oppSR.wave_type(r,:);

SRdirect=SourceReceiverPairs(myTransducers,myPlattens,submap,subwave);
fig_b=plotblockwithplattens(myBlock,myPlattens);

fig_handle=plotdirectrays(SRdirect,fig_b);

%% brute force automatic peakinf of first wave arrival using AIC
% Brice April 15 2019   -> to be  improved and then move to a function

% iterates on pair in the map
% we use the first acquisition sequence only here - Initial velocity. 
type=[];
slowness=zeros(SRdirect.n_pairs,3);
for p=1:SRdirect.n_pairs
    signal=squeeze(dataseq1(1,1:end,SRdirect.SRmap(p,1),SRdirect.SRmap(p,2)));
    
    type=[type; SRdirect.wave_type(p,:)];
    
    % figure
    % plot(signal,'.')
    % aic ?
    tic
    aic = AIC(signal);
    toc
    
    % figure
    % plot(aic)
    %
    % figure
    % plot(diff(aic))
    
    
    % now let's peak the arrival using 3 different criteria
    
    [f k]=min(aic(endnoise:end-10));  % minimum of AIC - not reliable (with multiple arrival)    
    [f2 k2]=max(diff(aic(endnoise:end-10))); % max of AIC derivatives (>0) -> good for first arrival in most cases    
    [f3 k3 ]=max(signal(endnoise:end)); % maximum signal peak -> lower velocity 

    dtpick=[k+endnoise k2+endnoise k3+endnoise];
    slowness(p,:)=dtpick.*dt/SRdirect.distances(p); % estimate slowness 1/Vel
    
end
vel=1./slowness; % corresponding velocity
% some plottings
figure
plot(vel);

% from the plot  - we can extract "manually" P and S waves 
% below is very specific to that case
dvp=[[1:30],[32:54],[56:length(vel(:,1))]];

dvs=[15,28,32,51,52,54,55];

vpm=mean(vel(dvp,2))
sqrt(var(vel(dvp,2)))

vsm=mean(vel(dvs,3))
sqrt(var(vel(dvs,3)))

gm=vsm^2*2700

km=vpm^2*2700-2./3*gm
ee=9*gm*km/(gm+3*km)
nu=(3*km-2*gm)/(2*(gm+3*km))

%% look at one particular combinaison of the source-receiver pair in the map
p=51;
disp(['We look at Source ',string(SRdirect.SRmap(p,1)),' and receiver ',string(SRdirect.SRmap(p,2)) ,' of type', string(SRdirect.wave_type(p,:)) , ' S-R distance (m):',string(SRdirect.distances(p))]);
signal=squeeze(dataseq1(1,1:end,SRdirect.SRmap(p,1),SRdirect.SRmap(p,2)));


figure 
title('Signal')
plot(signal,'.')
% aic ?
tic
aic = AIC(signal);
toc

figure
plot(aic)
%
figure
plot(diff(aic))

% Note - things to be done
%  picking S arrivals after the P ! -> how to automate etc.
%


%% velocity measurement via auto-correlation -> Thomas / Dmitry
% example goes here 

%% Trying to look at diffracted signals during the HF experiment

% first we chooise some candidates pairs 

% here we take all pairs that are neither intra platen neither opposite
% platen
oppSR = AllPairsOppositePlattens(myTransducers,myPlattens);
oppexceptIntra=AllexceptintraPairs(myTransducers,myPlattens);

[my_map r]=setdiff(oppexceptIntra.SRmap,oppSR.SRmap,'rows');


SRdiffAll=SourceReceiverPairs(myTransducers,myPlattens,my_map,oppexceptIntra.wave_type(r,:));
% takes the ones that a have snr above 2  on the first sequence (stringent)
snrv=zeros(length(SRdiffAll.SRmap),1);
for i=1:length(snrv)
    snrv(i)=normAll(SRdiffAll.SRmap(i,1),SRdiffAll.SRmap(i,2));
end


[r c ~]=find(snrv>2);

snrvf=snrv(r);
% create the new S-R objects 
submap= SRdiffAll.SRmap(r,:);
subwave=SRdiffAll.wave_type(r,:);
SRdiff=SourceReceiverPairs(myTransducers,myPlattens,submap,subwave);
% plot it
fig_b=plotblockwithplattens(myBlock,myPlattens);

fig_handle=plotdirectrays(SRdiff,fig_b);


%% compute integrated SNR for all sequences chosen for all S_R diffractor candidates pairs
endnoise = 200;
snrdt=zeros(SRdiff.n_pairs,1);
for p=1:SRdiff.n_pairs
    
    basesignal=squeeze(dataseq1(1,1:end,SRdiff.SRmap(p,1),SRdiff.SRmap(p,2)));
    allsignal=squeeze(dataseq1(1:end,1:end,SRdiff.SRmap(p,1),SRdiff.SRmap(p,2)));
    
    [dataSubs] = base_signal_substraction(allsignal,basesignal);
    snrdt(p)=sum(abs(dataSubs(:,endnoise:end).^2),'all');
    
end
figure 
plot(snrdt)

% here we choose the ones with snrdt > 20 (case dependendt) 
[kp]=find(snrdt>0.8)
%% Create the corresponding map 

submap= SRdiff.SRmap(kp,:);
subwave=SRdiff.wave_type(kp,:);

SR_i=SourceReceiverPairs(myTransducers,myPlattens,submap,subwave);
fig_b=plotblockwithplattens(myBlock,myPlattens);

fig_handle=plotdirectrays(SR_i,fig_b);

%% diffraction look at one candidate S_R pair
p= 2 
disp(["Diffraction"]);
disp(['We look at Source ',string(SR_i.SRmap(p,1)),' and receiver ',string(SR_i.SRmap(p,2)) ,' of type', string(SR_i.wave_type(p,:)) , ' S-R distance (m):',string(SR_i.distances(p))]);

basesignal=squeeze(dataseq1(1,endnoise:end,SR_i.SRmap(p,1),SR_i.SRmap(p,2)));
allsignal=squeeze(dataseq1(1:end,endnoise:end,SR_i.SRmap(p,1),SR_i.SRmap(p,2)));

[dataSubs] = base_signal_substraction(allsignal,basesignal);

figure;


wiggle(1e6*T(endnoise:end)',300*tt',dataSubs','+');  % here the x scale is not so nice

%wiggle( dataSubs','-');  % plot without scales

 
%% fluid thickness analysis - example of few pairs
% define parameters
myFluid = Fluid('glycerol');
% NOTE the API SHOULD CHANGE and use Isotropic Solid object or TISolid
% object ! 
mySolid.Vp = 6500;
mySolid.rho = 2750;
timeWindow = [38 44]*1E-6;
freqWindow = [0.5 1.1]*1E6;

% reload dataset 
dataseq2 = load_data(fbin,1:nq,[3 9 14],[3 9 14]);   %--- here this where the use of the SRmap in one SR pair object is useful -> to replace
ns2 = size(dataseq2,3);
nr2 = size(dataseq2,4);

% extract on-axis pairs only  
D = reshape(squeeze(dataseq2(:,:,:,:)),size(dataseq2,1),[],ns2*nr2); % flatten 3D array
Dd2 = squeeze(D(:,:,1:nr2+1:end)); % extract 'diagonal' for all sequences

% call for the fluid thickness minimization
w = fluid_thickness(Dd2,squeeze(Dd2(1,:,:)),T,mySolid,myFluid,...
    timeWindow,freqWindow);

%% plot analyzed width 
figure 
plot(AcqTime(1:200),w(1:200,:))

% needs to create meaningful plots - locate of frac width in space and time
% etc.



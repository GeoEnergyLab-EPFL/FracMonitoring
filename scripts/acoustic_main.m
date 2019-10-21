%% cleanup first, set global parameters
%close all
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

% % injection on slate
% datamonth = 03;
% dataday = 01;
% starttime = '121330';

% test on gabbro
datamonth = 03;
dataday = 14;
starttime = '093621';

% test transducers
% % NS platten
% datamonth = 03;
% dataday = 27;
% starttime = '114711';.
% % EW platten
% datamonth = 03;
% dataday = 27;
% starttime = '175106';
% % all platten
% datamonth = 04;
% dataday = 05;
% starttime = '141723';

% test aluminium
% datamonth = 01;
% dataday = 29;
% starttime = '133723';


% data folder name from experiment date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];

% extract header info from JSON file
fjson = [datapath datafold num2str(starttime) '.json'];
[jsonhdr,myTransducers,myPlattens,myBlock] = load_header(fjson);

%% load data from bin file
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

fbin = [datapath datafold num2str(starttime) '.bin'];
[dataseq1, nq] = load_data(fbin,1:1,1:32,1:32);

%% low rate data
[lowrateData, lowrateTime, lowrateHeader] = load_lowrate(fbin);

% plot injection flows
figure
plot(lowrateTime,lowrateData(:,11),lowrateTime,lowrateData(:,14))
axis tight
xlim([lowrateTime(1) lowrateTime(end)])
ylim([0 1])

% plot pressure
figure
plot((1:length(lowrateTime))/10,lowrateData(:,10)/1E3)%,(1:length(lowrateTime))/10,lowrateData(:,13)/1E3)
axis tight
xlabel('Time (s)')
ylabel('Pressure (MPa)')


%% plot all
figure
plot(T*1E6,squeeze(dataseq1(1,:,:)))
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('All transducer combinations')


%% plot transducer pairs
% longitudinal transducers
figure
hold on
for i_pair = [1:23 25:28] %min([ns nr])
    plot(T*1E6,squeeze(dataseq1(1,:,i_pair,i_pair)))
end
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('Longitudinal transducer pairs')

% shear transducers
figure
hold on
for i_pair = [24 29:32]
    plot(T*1E6,squeeze(dataseq1(1,:,i_pair,i_pair)))
end
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('Shear transducer pairs')

% plot arrivals for one source
i_seq = 1;
i_source = 4;
figure
imagesc(1:nr,T*1E6,squeeze(dataseq1(i_seq,:,i_source,:)))
xlabel('Receiver number')
ylabel('Time (\mus)')
title(['All receivers for source # ' num2str(i_source)])
caxis([-1 1]*0.01)

% plot arrivals for one receiver
i_seq = 1;
i_receiver = 24;
figure
imagesc(1:nr,T*1E6,squeeze(dataseq1(i_seq,:,:,i_receiver)))
xlabel('Receiver number')
ylabel('Time (\mus)')
title(['All sources for receiver # ' num2str(i_receiver)])
caxis([-1 1]*0.01)

%% look for good source-receiver pairs
D = reshape(squeeze(dataseq1(1,:,:,:)),[],ns*nr); % flatten 3D array
Dd = squeeze(D(:,1:nr+1:end)); % extract 'diagonal'
clearvars D

% remove excitation noise
endnoise = 100;
% L2 norm 
N = vecnorm(Dd(endnoise:end,:));
% plot L2 norm strength
figure
disp('plotting plot L2 norm strength')
bar(1:nr,N)
axis([0.5 nr+0.5 0 max(N)])
xlabel('Source-receiver pair')
ylabel('Signal strength')

%% all source receiver pair L2-norms
% extract L2 norms
Nall = squeeze(vecnorm(squeeze(dataseq1(1,endnoise:end,:,:)),2,1));
% plot
figure
imagesc(1:nr,1:ns,Nall)
axis square
caxis([0 1]*2.5)
colormap('jet')
xlabel('Receiver #')
ylabel('Source #')

hold on
plot(1:32,1:32,'.k')

%% plot for specific source or receiver
i_source = 28;
figure
imagesc(1:32,T*1E6,squeeze(dataseq1(1,:,i_source,:)))
xlabel('Transducer number')
ylabel('Time (\mus)')
caxis([-1 1]*0.02)

%% look into change in time for one pair
i_pair = 32;
datasource1 = load_data(fbin,1:nq,i_pair,i_pair);

figure
plot(T*1E6,squeeze(datasource1))

%% initial velocity analysis
dataset = Dd(:,1:16);
% reference signal read (without a sample), p-wave
dataRef = load_data([datapath '19-01-11/163814.bin'],1,1,1)';
Vp_test = velocity_measurements(dataset,dataRef,dt);


%% V ADVANCED start

%% Separating P and S from class, create a matrix  S | R | distance | type | angle

oppSR = AllPairsOppositePlattens(myTransducers,myPlattens); % change to not only opposite to not on the same platten !!! TO DO

VmatrixIn.SRmap = oppSR.SRmap; % S/R
VmatrixIn.Ttype = string(oppSR.wave_type); % type P/S
VmatrixIn.Length = oppSR.distances; % length
VmatrixIn.N_pairs = oppSR.n_pairs; % n_pairs
VmatrixIn.Angle = oppSR.directions; % angle??
VmatrixIn.dt = dt; %dt, tipestep b/w measurements

%% load data for V only useful SR pair /// CREATE FOR ALL SEQ TBD!!!
for iSR=1:VmatrixIn.N_pairs % for all SRmap
    VmatrixIn.dataset(:,iSR) = dataseq1(1,:,VmatrixIn.SRmap(iSR,1),VmatrixIn.SRmap(iSR,2)); % create a dataset in Vmatrix (dont comapre with D -- it is by receivers)
end

% Ref
% dataset = Dd(:,1:5); % only opposite pairs and first seq ???
% reference signal read (without a sample), p-wave
VmatrixIn.dataRef_P = load_data([datapath '19-01-11/163814.bin'],1,1,1)'; % Ref signal:Swave = 162209, Pwave = 163814
VmatrixIn.dataRef_S = load_data([datapath '19-01-11/162209.bin'],1,1,1)'; % Ref signal:Swave = 162209, Pwave = 163814

%% V velocity function

[VmatrixOut] = velocity_advanced(VmatrixIn, {30,50}, {60,80}, ["Full", "Part", "GI", "AIC"]); % Vmatrix, {Pwave Tmin, Pwave Tmax}, {S wave Tmin, Swave Tmax}, ["Full", "Part", "GI", "AIC", "STA_LTA"], dt


%% fluid thickness 
% define parameters
myFluid = Fluid('glycerol');
mySolid.Vp = 6500;
mySolid.rho = 2750;
timeWindow = [38 44]*1E-6;
freqWindow = [0.5 1.1]*1E6;

% extract on-axis pairs only
D = reshape(squeeze(dataseq1(:,:,:,:)),size(dataseq1,1),[],ns*nr); % flatten 3D array
Dd2 = squeeze(D(:,:,1:nr+1:end)); % extract 'diagonal' for all sequences

w = fluid_thickness_analysis(Dd2,squeeze(Dd2(1,:,:)),T,mySolid,myFluid,...
    timeWindow,freqWindow);


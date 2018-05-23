close all
clearvars
home

%% data path on ENACDrives
if isunix
    [~, uid] = system('id -u');
    [~, username] = system('whoami');
    datapath = ['/run/user/' uid(1:end-1) '/gvfs/smb-share:domain=INTRANET,server=enac1files.epfl.ch,'...
    'share=gel,user=' username(1:end-1) '/research/Experiments/HF-Experiments/GEL-data/'];
elseif ispc
    datapath = 'Z:/research/Experiments/HF-Experiments/GEL-data/';
end

datafold = '18-05-09/';
endfile = '134903';

%datafold = '18-04-19/';
%endfile = '211411';

%datafold = '18-04-26/';
%endfile = 132320;

% error if data folder not found
if ~exist([datapath datafold],'dir')
    error('acoustic_read:datafolder','The data folder does not exist')
end
% error if timing and params files not found
if ~exist([datapath datafold 'params_' num2str(endfile) '.txt'],'file')
    error('acoustic_read:params','The params file does not exist')
end
if ~exist([datapath datafold 'timing_' num2str(endfile) '.txt'],'file')
    error('acoustic_read:params','The timing file does not exist')
end

% extract timestamp list
fid = fopen([datapath datafold 'timing_' num2str(endfile) '.txt'],'r');
CellTimes = textscan(fid,'%s');
fclose(fid);
AcqTime = datetime(CellTimes{1},'Format','HH:mm:ss');  % time in HH:mm:ss
AcqTime.Year = 2018;
AcqTime.Month = 05;
AcqTime.Day = 08;
% filenum list
tmp = datevec(AcqTime);
filenum  = tmp(:,4)*10000+tmp(:,5)*100+tmp(:,6);

% get first acquisition time
[~, ~, ~, HH, mm, ss,] = datevec(datetime(CellTimes{1}{1}));
startfile = HH*1E4+mm*1E2+ss;

%% pressure profile
% load pressure data from two pressure gauges
fid = fopen([datapath datafold 'voltage_' num2str(endfile) '.txt'],'r');
CellTimes = textscan(fid,'%s %f %f');
fclose(fid);
PressureTime = datetime(CellTimes{1},'Format','HH:mm:ss');  % time in HH:mm:ss
Pressure = [60*CellTimes{2} 60*CellTimes{3}]; % pressure in bars

% plot both pressures in time
figure
disp('plotting pressure over time')
plot(PressureTime,Pressure)
%axis(datenum([PressureTime(1) PressureTime(end)]), [0 60]) %error
xlabel('Time')
ylabel('Pressure (bar)')

% difference between two gauges
figure
disp('plotting gauges difference')
plot(PressureTime,Pressure(:,2)-Pressure(:,1))
%axis([datenum([PressureTime(1) PressureTime(end)]) -1 1]) %error
xlabel('Time')
ylabel('Pressure (bar)')

clearvars CellTimes
                                   
%% variable definitions and data file listing
% get folder info and data file list from start and endfile info
FolderInfo = dir([datapath datafold]);
startindex = find(strcmp({FolderInfo.name}, ['data_' num2str(startfile,'%06d') '.bin']) == 1);
endindex = find(strcmp({FolderInfo.name}, ['data_' endfile '.bin']) == 1);

% read params file
% for now acquisition parameters are hard-coded, this will change with the
% updated acquisition VI from EPSLog and header info
ns = 8000;  % nb of time samples per acquisition
nt = 32;    % number of sources
nr = 32;    % number of receivers
nq = endindex-startindex+1;   % nb of acquisitions

Fs = 5E7;   % sampling frequency (hardware defined)
dt = 1/Fs;  % time step
t0 = 0;     % initial time
T = t0+dt*linspace(0,ns-1,ns)'; % time vector
Fn = 0.5*Fs;    % Nyquist frequency (Hz)

% DC filter definition
a = [1,-0.99];
b = [1,-1];

%% load initial acquisitions
% check the first nn acquisition sequences for a specific source
nn = 4;

% load data from bin files and DC filter it too
datatmp = zeros(nn,ns,nr*nt);
datafilt = zeros(size(datatmp));
for ii = 1:nn
    fid = fopen([datapath datafold FolderInfo(startindex+ii-1).name],'r');
    datatmp(ii,:,:) = fread(fid,[ns,nt*nr],'double');
    fclose(fid);
    datafilt(ii,:,:) = filtfilt(b,a,squeeze(datatmp(ii,:,:)));
end

% reshape data
dataInit2 = reshape(datatmp,nn,ns,nr,nt);
dataInit3 = reshape(datafilt,nn,ns,nr,nt);
% 1st index, nn, is nb of initial acquisitions to look at
% 2nd index, ns, is nb of time samples
% 3rd index, nr, in nb of receivers
% 4th index, nt, is nb of sources

clearvars datatmp datafilt

%% plot them
% time plot
jj = 4;
figure
disp('plotting source-receiver amplitude over time')
plot(T*1E6,dataInit2(:,:,jj,jj),T*1E6,dataInit3(:,:,jj,jj))
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
title(['source-receiver #' num2str(jj)])

% image plot
kk = 4;
figure, imagesc(0:nr-1,T*1E6,squeeze(dataInit3(1,:,:,kk)))
disp('plotting source #4 over time')
caxis([-1 1]*0.002)
colormap('jet')
colorbar
axis([0 nr-1 0 150])
xlabel('Receiver number')
ylabel('Time (\mus)')
title(['Source number ' num2str(kk)])

%% look for good source-receiver pairs
%D = bsxfun(@times, eye(size(data0(:,:,1))), data0);
D = reshape(squeeze(dataInit2(1,:,:,:)),[],nr*nt); % flatten 3D array
Dd = squeeze(D(:,1:nr+1:end)); % extract 'diagonal'
clearvars D

% check all on-axis pairs
figure
disp('plotting all on-axis pairs')
plot(T*1E6,Dd)

% remove excitation noise
endnoise = 100;
% L2 norm 
N = sqrt(sum(Dd(endnoise:end,:).^2,1));
% plot L2 norm strength
figure
disp('plotting plot L2 norm strength')
bar(1:nr,N)
axis([0.5 nr+0.5 0 1])
xlabel('Source-receiver pair')
ylabel('Signal strength')

%% save selected pairs in new array
% Pairs = {3,11,17,23};
% DataS = zeros(nq,ns,length(Pairs));
% for ii = 1:length(Pairs)
%     DataS(:,:,ii) = squeeze(data3(:,:,Pairs{ii},Pairs{ii}));
% end

%% clear initial data
clearvars dataInit2 dataInit3

%% load selected source receiver pair for all acquisition sequences
% select pair
jj = 14; % receiver number
kk = 14; % source number
% load data from bin files and DC filter it too
dataPair = zeros(ns,nq);
dataPairFilt = zeros(size(dataPair));
for ii = 1:nq
    fid = fopen([datapath datafold FolderInfo(startindex+ii-1).name],'r');
    fseek(fid,(kk-1)*ns*nr*8+((jj-1)*ns*8),'bof'); 
    dataPair(:,ii) = fread(fid,ns,'double');
    fclose(fid);
    dataPairFilt(:,ii) = filtfilt(b,a,squeeze(dataPair(:,ii)));
end

%% analyze waveform changes from fluid layer
% quick figure for selected source-receiver pair
figure
disp('plotting figure for selected source-receiver pair')
plot(T*1E6,dataPair)
axis([80 120 [-1 1]*5E-1])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
%celllgd = cellstr(num2str(filenum','%4d'));
%legend(celllgd{:})

% lowpass filter
flowpass = 4E6; % cut at 4 MHz
[b, a] = butter(2,flowpass/Fn);
dataLowFilt = filtfilt(b,a,dataPair);

% quick figure
kk = 20;
figure
disp('plotting next figure for selected source-receiver pair')
plot(T*1E6,dataLowFilt(:,1:kk))
axis([80 120 [-1 1]*1])
%axis([88 112 [-1 1]*5E-1])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
celllgd = cellstr(AcqTime);
legend(celllgd{1:kk})

% max as function of time
ampmax = max(abs(hilbert(dataLowFilt)));
figure, plot(AcqTime, ampmax,'o-')
disp('plotting max amplitude of time')
%axis([datenum([AcqTime(1) AcqTime(end)]) [0 8]*1E-1])
xlabel('Time')
ylabel('Max Amplitude')

%% freq domain for thickness estimation
% move to frequency domain
nfft = 2^nextpow2(ns);
U = fft(dataPair,nfft);
U = U(1:nfft/2+1,:)/nfft;
freq = Fn*linspace(0,1,nfft/2+1)';

% graphical check
figure
disp('plotting frequency vs amplitude')
plot(freq*1E-6,abs(U))
xlabel('Frequency (MHz)')
ylabel('Amplitude (a.u.)')
axis([0 2 0 5E-3])

%% thickness estimation
% experimental constants
rho_plexi = 1180;   % density PMMA
v_plexi = 2790;     % velocity PMMA
rho_glycerol = 1260;% density glycerol
v_glycerol = 1960;  % velocity glycerol
% choose the right solid and fluid properties
solid = 'plexi';
fluid = 'glycerol';
% set the fluid and solid properties for the selected materials
rho_solid = eval(['rho_' solid]);
v_solid = eval(['v_' solid]);
rho_fluid = eval(['rho_' fluid]);
v_fluid = eval(['v_' fluid]);

% frequency band to integrate over
flow = find(freq>=.6E6,1);
fhigh = find(freq>=.9E6,1);
freqband = freq(flow:fhigh);

% transmission coef
h = (-250:0.5:250)*1E-6;
alpha = 2*pi*freqband*h/v_fluid;  % freq * thickness
Zr = rho_fluid*v_fluid/(rho_solid*v_solid);
rff = (Zr-1)/(Zr+1);
Trans = ((1-rff^2)*exp(-1i*alpha))./(1-rff^2*exp(-2*1i*alpha));

% objective function
Fun = zeros(nq,length(h));
for ii = 1:nq
    Fun(ii,:) = sum(abs((U(flow:fhigh,ii)*ones(size(h))-(U(flow:fhigh,1)...
        *ones(size(h)).*Trans)).^2),1);
end

% locate min
[~, hmin] = min(Fun,[],2);
figure, plot(AcqTime(1:end),h(hmin)*1E6,'*:')
disp('plotting fluid layer thikness vs time')
xlabel('Time')
ylabel('Fluid layer thickness (\mum)')
%axis([datenum([AcqTime(1) AcqTime(end)]) [-1 1]*200])

% % plots
% figure
% plot(h*1E6,Fun(end,:))
% xlabel('Thickness (\mum)')
% ylabel('Objective function')

%% thickness estimation for all raypaths normal to the fracture plane (vertical)
% define source-receiver paris to consider
SRpairs = 1:16;
% initialize objective function
Fun = zeros(length(SRpairs),nq,length(h));
% loop on source-receiver pairs
for jj = 1:length(SRpairs)
    % get freq data from first acquisition for reference
    fid = fopen([datapath datafold FolderInfo(startindex).name],'r');
    fseek(fid,(SRpairs(jj)-1)*ns*nr*8+((SRpairs(jj)-1)*ns*8),'bof');
    datatmp = fread(fid,ns,'double');
    fclose(fid);
    datatmpFilt = filtfilt(b,a,datatmp);
    Utmp = fft(datatmp,nfft);
    Uref = Utmp(1:nfft/2+1)/nfft;
    % now loop on remaining acquisitions to compare
    for ii = 2:nq
        fid = fopen([datapath datafold FolderInfo(startindex+ii-1).name],'r');
        fseek(fid,(SRpairs(jj)-1)*ns*nr*8+((SRpairs(jj)-1)*ns*8),'bof');
        datatmp = fread(fid,ns,'double');
        fclose(fid);
        datatmpFilt = filtfilt(b,a,datatmp);
        Utmp = fft(datatmp,nfft);
        U = Utmp(1:nfft/2+1)/nfft;
        Fun(jj,ii,:) = sum(abs((U(flow:fhigh)*ones(size(h))-...
            (Uref(flow:fhigh)*ones(size(h)).*Trans)).^2),1);
    end
end

%% locate and plot estimated thicknesses from min of objective function
[~, hmin] = min(Fun,[],3);
figure, plot(AcqTime(1:end),h(hmin)*1E6)
disp('plotting fluid layer thikness vs time (2)')
xlabel('Time')
ylabel('Fluid layer thickness (\mum)')
%axis([datenum([AcqTime(1) AcqTime(end)]) [-2.5 0.5]*100])
legend(cellstr(num2str(SRpairs','%d')));

%% add events to plot
% define times
EventTimes = datetime({'10:00:00','10:10:00','10:14:00','10:30:00',...
    '11:30:00','11:54:00'},'Format','HH:mm:ss');
EventTimes.Year = 2018;
EventTimes.Month = 05;
EventTimes.Day = 08;
% define vertical range
VertRange = [-1 1]*1E3;
hold on
for ii = 1:length(EventTimes)
    plot([EventTimes(ii) EventTimes(ii)],VertRange,':k')
end

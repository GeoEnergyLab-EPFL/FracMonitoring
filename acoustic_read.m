close all
clearvars
home

debugmode = 'off';

%% data path on ENACDrives
if isunix
    [~, uid] = system('id -u');
    [~, username] = system('whoami');
    datapath = ['/run/user/' uid(1:end-1) '/gvfs/smb-share:domain=INTRANET,server=enac1files.epfl.ch,'...
    'share=gel,user=' username(1:end-1) '/research/Experiments/HF-Experiments/GEL-data/'];
elseif ispc
    datapath = 'Y:/research/Experiments/HF-Experiments/GEL-data/';
end

datafold = '18-04-19/';
endfile = '211411';

%datafold = '18-04-26/';
%endfile = 132320;

% error if data folder not found
if ~exist([datapath datafold],'dir')
    error('acoustic_read:datafolder','The data folder doesn''t exist')
end
% error if timing and params files not found
if ~exist([datapath datafold 'params_' num2str(endfile) '.txt'],'file')
    error('acoustic_read:params','The params file doesn''t exist')
end
if ~exist([datapath datafold 'timing_' num2str(endfile) '.txt'],'file')
    error('acoustic_read:params','The timing file doesn''t exist')
end

% extract timestamp list
fid = fopen([datapath datafold 'timing_' num2str(endfile) '.txt'],'r');
C = textscan(fid,'%s');
fclose(fid);
AcqTime = datetime(C{1},'Format','HH:mm:ss');  % time in HH:mm:ss
% filenum list
tmp = datevec(AcqTime);
filenum  = tmp(:,4)*10000+tmp(:,5)*100+tmp(:,6);

% get first acquisition time
[~, ~, ~, HH, mm, ss,] = datevec(datetime(C{1}{1}));
startfile = HH*1E4+mm*1E2+ss;
                                   
%% variable definitions and data file listing
% get folder info and data file list from start and endfile info
FolderInfo = dir([datapath datafold]);
startindex = find(strcmp({FolderInfo.name}, ['data_' num2str(startfile) '.bin']) == 1);
endindex = find(strcmp({FolderInfo.name}, ['data_' num2str(endfile) '.bin']) == 1);

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
plot(T*1E6,dataInit2(:,:,jj,jj),T*1E6,dataInit3(:,:,jj,jj))
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
title(['source-receiver #' num2str(jj)])

% image plot
kk = 4;
figure, imagesc(0:nr-1,T*1E6,squeeze(dataInit3(1,:,:,kk)))
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

% check all on-axis pairs
figure
plot(T*1E6,Dd)

% remove excitation noise
endnoise = 100;
% L2 norm 
N = sqrt(sum(Dd(endnoise:end,:).^2,1));
% plot L2 norm strength
figure
bar(1:nr,N)
axis([0.5 nr+0.5 0 1])
xlabel('Source-receiver pair')
ylabel('Signal strength')

%% save selected pairs in new array
Pairs = {3,11,17,23};
DataS = zeros(nq,ns,length(Pairs));
for ii = 1:length(Pairs)
    DataS(:,:,ii) = squeeze(data3(:,:,Pairs{ii},Pairs{ii}));
end

%% pressure profile
fid = fopen([datapath datafold 'voltage_' num2str(endfile) '.txt'],'r');
C = textscan(fid,'%s %f %f');
fclose(fid);
PressureTime = datetime(C{1},'Format','HH:mm:ss');  % time in HH:mm:ss
Pressure = [60*C{2} 60*C{3}]; % pressure in bars

figure
plot(PressureTime,Pressure)
axis([datenum([PressureTime(1) PressureTime(end)]) 0 220])
xlabel('Time')
ylabel('Pressure (bar)')

% difference between two gauges
figure
plot(PressureTime,Pressure(:,2)-Pressure(:,1))
axis([datenum([PressureTime(1) PressureTime(end)]) 0 1])
xlabel('Time')
ylabel('Pressure (bar)')

%% load selected source receiver pair for all acquisition sequences
% select pair
jj = 3;
% load data from bin files and DC filter it too
dataPair = zeros(ns,nq);
dataPairFilt = zeros(size(datatmp));
for ii = 1:nq
    fid = fopen([datapath datafold FolderInfo(startindex+ii-1).name],'r');
    fseek(fid,(jj-1)*ns*nr*8+((jj-1)*ns*8),'bof'); 
    dataPair(:,ii) = fread(fid,ns,'double');
    fclose(fid);
    dataPairFilt(:,ii) = filtfilt(b,a,squeeze(dataPair(:,ii)));
end
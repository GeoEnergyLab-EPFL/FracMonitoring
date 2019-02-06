close all
clearvars
home

%% data path on ENACDrives
%datapath = pathbyarchitecture;
datapath = '/home/tblum/data/';

% 2018 acquisitions
datayear = 18;

datamonth = 09; dataday = 18; endtime = '121809'; % injection on marble block 2
%datamonth = 09; dataday = 18; endtime = '100103'; % test before injection on marble block 2

%datamonth = 09; dataday = 14; endtime = '190223'; % injection on marble block 2
%datamonth = 09; dataday = 14; endtime = '123341'; % test before injection on marble block 2

%datamonth = 08; dataday = 21; endtime = '190223'; % injection on marble block 2
%datamonth = 08; dataday = 21; endtime = '121046'; % test before injection on marble block 2

%datamonth = 08; dataday = 17; endtime = '180924'; % fracture in marble block
%datamonth = 08; dataday = 17; endtime = '102644'; % Marble block with amp, sec test before injection
%datamonth = 08; dataday = 16; endtime = '174852'; % Marble block with amp, test before injection

%datamonth = 07; dataday = 26; endtime = '134021'; % Slate block with pulser, test before injection

%datamonth = 07; dataday = 24; endtime = '155252'; % Slate block with pulser, inside cell with ground connection (T-B only)
%datamonth = 07; dataday = 24; endtime = '131632'; % Slate block with pulser, inside cell but isolated (T-B only)
%datamonth = 07; dataday = 24; endtime = '111240'; % Slate block with pulser, taken out of the cell (T-B only)
%datamonth = 07; dataday = 24; endtime = '091201'; % Slate block with pulser, inside the cell (T-B only)
%datamonth = 07; dataday = 23; endtime = '173517'; % Slate block with pulser, outside the cell (T-B only)

%datamonth = 07; dataday = 23; endtime = '121117'; % Slate block with pulser again, no fracture

%datamonth = 07; dataday = 18; endtime = '163738'; % Slate block with pulser again, no fracture
%datamonth = 07; dataday = 17; endtime = '131713'; % Slate block with pulser

%datamonth = 07; dataday = 09; endtime = '173250'; % PMMA half block with pulser for debugging
%datamonth = 07; dataday = 04; endtime = '161650'; % acquisition on PMMA half block for debugging
%datamonth = 07; dataday = 04; endtime = '174557'; % acquisition on PMMA half block with 16 ch only

%datamonth = 06; dataday = 05; endtime = '175813'; % test acquisition cement block C

%datamonth = 05; dataday = 28; endtime = '170423'; % 

%datamonth = 05; dataday = 28; endtime = '100142'; % 

%datamonth = 05; dataday = 08; endtime = '175209'; % test PMMA half blocks

%datamonth = 05; dataday = 09; endtime = '134903'; % injection PMMA half blocks

%datamonth = 04; dataday = 19; endtime = '211411'; % injection in cement block A

%datamonth = 03; dataday = 26; endtime = '163010'; % acquisition test on aluminium block

% 2017 acquisitions
%datayear = 17;

%datamonth = 11; dataday = 09; endtime = '121428'; % static test on single PMMA half block, ONLY 16 CHNLS

% build datafold from date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];
% error if data folder not found
if ~exist([datapath datafold],'dir')
    error('acoustic_read:datafolder','The data folder does not exist')
end
% error if timing and params files not found
if ~exist([datapath datafold 'params_' num2str(endtime) '.txt'],'file')
    error('acoustic_read:params','The params file does not exist')
end
if ~exist([datapath datafold 'timing_' num2str(endtime) '.txt'],'file')
    error('acoustic_read:params','The timing file does not exist')
end

% extract timestamp list
fid = fopen([datapath datafold 'timing_' num2str(endtime) '.txt'],'r');
CellTimes = textscan(fid,'%s');
fclose(fid);
AcqTime = timestamps(CellTimes{1},datayear,datamonth,dataday);

% filenum list
tmp = datevec(AcqTime);
filenum  = tmp(:,4)*10000+tmp(:,5)*100+tmp(:,6);

% get first acquisition time
startfile = tmp(1,4)*1E4+tmp(1,5)*1E2+tmp(1,6);

clearvars tmp

%% pressure profile
% load pressure data from two pressure gauges
fid = fopen([datapath datafold 'voltage_' num2str(endtime) '.txt'],'r');
CellTimes = textscan(fid,'%s %f %f');
fclose(fid);
PressureTime = timestamps(CellTimes{1},datayear,datamonth,dataday);
Pressure = 6000*[CellTimes{2} CellTimes{3}]; % pressure in MPa

% plot both pressures in time
figure
disp('plotting pressure over time')
plot(PressureTime,Pressure*1E-3) % change to MPa
xlim([PressureTime(1) PressureTime(end)])
ylim([0 40])
xlabel('Time')
ylabel('Pressure (MPa)')

% % difference between two gauges
% figure
% disp('plotting gauges difference')
% plot(PressureTime,Pressure(:,2)-Pressure(:,1))
% xlim([PressureTime(1) PressureTime(end)])
% ylim([-1 1]*100)
% xlabel('Time')
% ylabel('Pressure (kPa)')

clearvars CellTimes

%% Look at pressure changes and injection flow
% basic derivation
figure
plot(PressureTime(1:end-1),diff(Pressure))
xlim([PressureTime(500) PressureTime(end)])
ylim([-1 1]*50)
xlabel('Time')
ylabel('dP/dt (kPa/s)')

% filter first
fnpressure = 0.5;   % Nyquist freq for pressure data (0.5 Hz)
forder = 2; % filter order
flow = 0.2; % low-pass freq
[B, A] = butter(forder,flow/fnpressure,'low');
PressureFilt = filtfilt(B,A,Pressure);

% other filter
PressureMed = medfilt1(Pressure);

% % remove spikes and then filter
% thres = 100; % threshold pressure noise jump in kPa
% Idx = find(diff(Pressure(:,1))>thres);  % first find positive jump
% Idx = Idx((Pressure(Idx+1,1)-Pressure(Idx,1))>thres); % then negative jump behind
% % set spikes to NaN
% Pressure2 = Pressure;
% Pressure2(Idx+1) = NaN;
% % then interpolate
% Pressure2 = interp1(Pressure(~(Idx-1)),Pressure2(~(Idx-1)),(Idx-1)');
% Pressure2 = medfilt1(Pressure2);

% try other method
Pressure3 = filloutliers(Pressure,'center','movmedian',5,1);

% plot filtered pressures in time
figure
plot(PressureTime,Pressure3*1E-3) % change to MPa
xlim([PressureTime(1) PressureTime(end)])
ylim([0 40])
xlabel('Time')
ylabel('Pressure (MPa)')

% improved derivation
figure
plot(PressureTime(1:end-1),diff(Pressure3))
xlim([PressureTime(500) PressureTime(end)])
ylim([-1 1]*50)
xlabel('Time')
ylabel('dP/dt (kPa/s)')

%% make nice fig for poster
% find time of first flow
hh = 10; %11;
mm = 20; %42;
t1flow = find(PressureTime.Hour>=hh&PressureTime.Minute>=mm,1);
% time of lower flow rate
mm = 48;
t2flow = find(PressureTime.Hour>=hh&PressureTime.Minute>=mm,1);
% time of initiation
hh = 12; %14;
mm = 00; %19;
tinit = find(PressureTime.Hour>=hh&PressureTime.Minute>=mm,1);
% time of atm pressure
hh = 12;
mm = 12;
tatm = find(PressureTime.Hour>=hh&PressureTime.Minute>=mm,1);

% make figure
figure
set(gcf,'Position',[2286,1,800,400])
plot(PressureTime,Pressure3*1E-3) % change to MPa
% plot vertical markers
hold on
plot(PressureTime([t1flow t1flow]),[0 40],'k:')
%plot(PressureTime([t2flow t2flow]),[0 40],'k:')
plot(PressureTime([tinit tinit]),[0 40],'k:')
plot(PressureTime([tatm tatm]),[0 40],'k:')

datetick('x',15)
xlim([PressureTime(1) PressureTime(end)])
ylim([0 40])
% annotations
xlabel('Time')
ylabel('Pressure (MPa)')
title('Pressurization curve')
%text(PressureTime(1),5,'Preparation')
text(PressureTime(1),35,'High flow')
text(PressureTime(2000),35,'Lower injection flow')
text(PressureTime(6000),25,{'Breakdown and','depressurization'})

%% log-log pressure decrease
hh = 18;
mm = 09;

%hh = 14;
%mm = 19;
datetime0loglog = PressureTime(1)+hours(hh-PressureTime(1).Hour)+minutes(mm-PressureTime(1).Minute);
t0loglog = find(PressureTime>=datetime0loglog,1);

figure
loglog(0:length(Pressure(t0loglog:end,1))-1,Pressure(t0loglog:end,1)*1E3)

figure
plot(log10(0:length(Pressure(t0loglog:end,1))-1),log10(Pressure(t0loglog:end,1)*1E3)) % change to MPa
xlabel('log(Time)')

%% variable definitions and data file listing
% get folder info and data file list from start and endfile info
FolderInfo = dir([datapath datafold]);
startindex = find(strcmp({FolderInfo.name}, ['data_' num2str(startfile,'%06d') '.bin']) == 1);
endindex = find(strcmp({FolderInfo.name}, ['data_' endtime '.bin']) == 1);

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

%% load selected acquisitions
% check the first qq acquisition sequences after q0
q0 = 2;
qq = 2;
% % if 16 ch (for 2017 only)
% nt = 16;
% nr = 16;

[dataInit2, dataInit3] = loadbindata(ns,nr,nt,nq,q0,qq,[datapath datafold],FolderInfo,startindex,a,b);

%% plot them
% time plot
jj = 2; % source-receiver pair
figure
disp(['plotting source-receiver #' num2str(jj) ' amplitude over time'])
plot(T*1E6,dataInit2(:,:,jj,jj),T*1E6,dataInit3(:,:,jj,jj))
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
title(['source-receiver #' num2str(jj)])

% image plot
kk = 7; % source number
figure, imagesc(0:nr-1,T*1E6,squeeze(dataInit3(1,:,:,kk)))
disp(['plotting receivers for source #' num2str(kk) 'over time'])
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
axis([0.5 nr+0.5 0 max(N)])
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
jj = 18; % receiver number
kk = 16; %24; % source number
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

%% plot of all signals for this pair
figure
plot(T*1E6,dataPair(:,2:end))
axis([[T(1) T(end)]*1E6 [-1 1]*0.4])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')

%% repeatability bug analysis
% XC = zeros(ns*2-1,nq);
% for ii = 1:nq
%     XC(:,ii) = xcorr(dataPair(:,2),dataPair(:,ii));
% end
% MX = max(XC,[],1);
% 
% figure
% plot(1:nq,MX)
% 
% IDX = find(MX>30);

%% look at changes from diffraction
% plot traces with time
figure
imagesc(1:nq,T*1E6,dataPair)
caxis([-1 1]*0.002)
%colormap('jet')
colorbar
xlim([1000 nq])
ylim([T(1) T(3000)]*1E6)
xlabel('Acquisition number')
ylabel('Time (\mus)')

% compute difference
navg = 800;
dataAvg = mean(dataPair(:,1:navg),2);
dataDiff = bsxfun(@minus,dataPair,dataAvg);
dataAvgFilt = mean(dataPairFilt(:,1:navg),2);
dataDiffFilt = bsxfun(@minus,dataPairFilt,dataAvgFilt);

% plot difference
difffig = figure;
imagesc(1:nq,T*1E6,dataDiffFilt)
%colormap('gray')
caxis([-1 1]*0.0005)
xlim([950 nq])
ylim([T(1000) T(2500)]*1E6)
xlabel('Acquisition number')
ylabel('Time (\mus)')
%print(difffig,['./figures/diffraction' num2str(kk) '-' num2str(jj)],'-dpng')

%% automated diffraction picking (for pair 7-24)
% first additional filtering
Flow = 1.5E6;   % lowpass freq
[bbutter, abutter] = butter(2,Flow/Fn,'low');
dataButter = filtfilt(bbutter,abutter,dataPairFilt);
dataAvgButter = mean(dataButter(:,1:navg),2);
dataDiffButter = bsxfun(@minus,dataButter,dataAvgButter);
% plot difference
fbutter = figure;
imagesc(1:nq,T*1E6,dataDiffButter)
caxis([-1 1]*0.0005)
xlim([1000 nq])
ylim([T(1) T(3000)]*1E6)
xlabel('Acquisition number')
ylabel('Time (\mus)')

% guessed time from manual pick
q_guess = 1030; % trace number for first guess
t_guess = find(T>=40E-6,1); % guessed time in micros

% plot
figure
plot(T*1E6,dataDiff(:,q_guess),T*1E6,dataDiffFilt(:,q_guess),T*1E6,dataDiffButter(:,q_guess))
xlim([0 60])

% simple neighbor cross-correlation
normwin = find(T>=2.5E-6,1); % window half-length is 2 microsec
% improve initial guess by looking for peak
[~, t_tmp] = max(abs(hilbert(dataDiffButter(t_guess-normwin:t_guess+normwin,q_guess))));
t_inittrace = t_tmp+t_guess-normwin-1;
hold on, plot(T*1E6,abs(hilbert(dataDiffButter(:,q_guess))),'k')
tmptrc = abs(hilbert(dataDiffButter(:,q_guess)));
plot(T(t_inittrace)*1E6,tmptrc(t_inittrace),'*k')

% 
xcshift = zeros(1,nq);
for ii = q_guess+1:nq
    xc1 = dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii-1)/...
        max(abs(dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii-1)));
    xc2 = dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii)/...
        max(abs(dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii)));
    [xcres, lags] = xcorr(xc1,xc2);
    [~, yy] = max(xcres);
    xcshift(ii) = xcshift(ii-1)-(yy-1+lags(1));
end

figure(fbutter)
hold on
plot(1:nq,T(xcshift+t_inittrace)*1E6,'.k')


% figure of signal envelope
figure
imagesc(1:nq,T*1E6,abs(hilbert(dataDiffButter)))
caxis([0 1]*0.002)
xlim([1000 nq])
ylim([T(1) T(3000)]*1E6)
xlabel('Acquisition number')
ylabel('Time (\mus)')

%% peak detection from envelope and xcorr
% find max from envelope
tstart = find(T>=25E-6,1);
[pks,locs] = findpeaks(abs(hilbert(dataDiffButter(:,q_guess))));
[~, idx] = min(abs(locs-t_guess));
locstart = locs(idx);

% init pksshift
pksshift = locstart*ones(1,nq);
for ii = q_guess+1:nq
    [pks,locs] = findpeaks(abs(hilbert(dataDiffButter(:,ii))));
    [~, idxtmp] = min(abs(locs-pksshift(ii-1))); % find nearest peak
    pksshift(ii) = locs(idxtmp);
    if locs(idxtmp)>=pksshift(ii-1)+25&&idxtmp>1
        pksshift(ii) = locs(idxtmp-1);
    end
end

hold on
plot(1:nq,T(pksshift)*1E6,'w')

% pb after 1146
% guessed time from manual pick
q_guess = 1146; % trace number for first guess
t_guess = pksshift(q_guess); % guessed time in micros

% plot
figure
plot(T*1E6,dataDiff(:,q_guess),T*1E6,dataDiffFilt(:,q_guess),T*1E6,dataDiffButter(:,q_guess))
xlim([0 60])

% simple neighbor cross-correlation
normwin = find(T>=2.5E-6,1); % window half-length is 2 microsec
% improve initial guess by looking for peak
[~, t_tmp] = max(abs(hilbert(dataDiffButter(t_guess-normwin:t_guess+normwin,q_guess))));
t_inittrace = t_tmp+t_guess-normwin-1;
hold on, plot(T*1E6,abs(hilbert(dataDiffButter(:,q_guess))),'k')
tmptrc = abs(hilbert(dataDiffButter(:,q_guess)));
plot(T(t_inittrace)*1E6,tmptrc(t_inittrace),'*k')

% 
xcshift = zeros(1,nq);
for ii = q_guess+1:nq
    xc1 = dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii-1)/...
        max(abs(dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii-1)));
    xc2 = dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii)/...
        max(abs(dataDiffButter(t_inittrace+xcshift(ii-1)-normwin:t_inittrace+xcshift(ii-1)+normwin,ii)));
    [xcres, lags] = xcorr(xc1,xc2);
    [~, yy] = max(xcres);
    xcshift(ii) = xcshift(ii-1)-(yy-1+lags(1));
end

figure(fbutter)
hold on
plot(1:nq,T(xcshift+t_inittrace)*1E6,'.k')

%% Test directly from first max
% guessed time from manual pick
q_guess = 1026; % trace number for first guess
t_guess = find(T>=40E-6,1); % guessed time in micros
% t_guess = find(T>=39.2E-6,1); % for pair 7-24

% simple neighbor cross-correlation
frontwin = find(T>=1E-6,1);
backwin = find(T>=0.5E-6,1);
% improve initial guess by looking for peak
[~, t_tmp] = max(dataDiffButter(t_guess-frontwin:t_guess+backwin,q_guess));
t_inittrace = t_tmp+t_guess-frontwin-1;

% 
maxshift = zeros(1,nq);
ii1 = 1049; % forced point
ii2 = 1163; % forced point
ii3 = 1172;
for ii = q_guess+1:nq
    [~, idxtmp] = max(dataDiffButter(t_inittrace+maxshift(ii-1)-frontwin:...
        t_inittrace+maxshift(ii-1)+backwin,ii));
    maxshift(ii) = maxshift(ii-1)+idxtmp-frontwin-1;
%     if ii==ii1
%         maxshift(ii) = find(T>=41.3E-6,1)-t_inittrace-1;
%         %maxshift(ii) = 1726-t_inittrace-1; % for 7-24
%     end
%     if ii==ii2
%         maxshift(ii) = find(T>=30.64E-6,1)-t_inittrace-1;
%         % maxshift(ii) = 1486-t_inittrace-1; % for 7-24
%     end
%     if ii==ii3
%        maxshift(ii) = find(T>=29.82E-6,1)-t_inittrace-1;
%     end
end

%% make nice figure for poster
figure
set(gcf,'Position',[2286,1,1000,400])
imagesc(1:nq,T*1E6,dataDiffButter)
caxis([-1 1]*0.0005)
xlim([950 nq])
ylim([T(999) T(2500)]*1E6)
xlabel('Acquisition time')
ylabel('Travelime (\mus)')
xticklabels(cellstr(AcqTime(xticks)))
hold on
plot(1:nq,T(maxshift+t_inittrace)*1E6,'r')
title(['Diffraction arrival pair ' num2str(kk) '-' num2str(jj)'])

% add forced points
% plot(ii1,T(maxshift(ii1)+t_inittrace)*1E6,'ro')
% plot(ii2,T(maxshift(ii2)+t_inittrace)*1E6,'ro')

%% save data to file
Datasave = [(1000:nq)', T(maxshift(1000:nq)+t_inittrace)*1E6];
csvwrite(['~/matlab/EPFL/diffraction_18-09-18/' num2str(kk) '-' num2str(jj)' '.csv'],...
    Datasave);

    
%% inversion diffraction 
x = 0:1:125;
halflength = 125;
fracdepth = 20;
sourcepos = 54;
receiverpos = 50;
dist1 = sqrt((halflength+fracdepth)^2+(sourcepos-x).^2);
dist2 = sqrt((halflength-x).^2+(receiverpos-fracdepth)^2);
v_marble = 6600;
traveltime = (dist1+dist2)*1E-3/v_marble;

vq = interp1(traveltime,x,T(maxshift+t_inittrace));
vq(isnan(vq)) = max(vq);
figure
plot(AcqTime(1000:end),vq(1000:end))
datetick('x',15)

%% add to thickness figure
yyaxis right
plot(AcqTime(1000:end),vq(1000:end))
ylabel('Fracture tip position (mm)')

%% analyze waveform changes from fluid layer
% quick figure for selected source-receiver pair
figure
disp('plotting figure for selected source-receiver pair')
plot(T*1E6,dataPair)
axis([50 70 [-1 1]*5E-1])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
%celllgd = cellstr(num2str(filenum','%4d'));
%legend(celllgd{:})

% lowpass filter
flowpass = 4E6; % cut at 4 MHz
[b, a] = butter(2,flowpass/Fn);
dataLowFilt = filtfilt(b,a,dataPair);

% quick figure
ii = 1;
figure
disp('plotting next figure for selected source-receiver pair')
plot(T*1E6,dataLowFilt(:,ii))
axis([30 50 [-1 1]*1E-1])
%axis([88 112 [-1 1]*5E-1])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
celllgd = cellstr(AcqTime);
legend(celllgd{ii})

% max as function of time
ampmax = max(abs(hilbert(dataLowFilt)));
figure
disp('plotting max amplitude of time')
plot(AcqTime, ampmax,'o-')
xlim([PressureTime(1) PressureTime(end)])
ylim([0 1]*3E-1)
xlabel('Acquisition time')
ylabel('Max Amplitude (a.u.)')

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
axis([0 2 0 2E-3])

%% thickness estimation
% experimental constants
rho_plexi = 1180;   % density PMMA
v_plexi = 2790;     % velocity PMMA
rho_glycerol = 1260;% density glycerol (from paper)
v_glycerol = 1960;  % velocity glycerol (from paper)
rho_cement = 2000;  % density cement (from measurement with Lionel Sofia)
v_cement = 4600;    % velocity cement (crude estimate from first arrival)
rho_marble = 2700;  % density marble (from literature)
v_marble = 6600;    % velocity marble (crude estimate from first arrival)
rho_silicone = 1050;% 
v_silicone = 1350;  %
% choose the right solid and fluid properties
solid = 'marble';
fluid = 'silicone';
% set the fluid and solid properties for the selected materials
rho_solid = eval(['rho_' solid]);
v_solid = eval(['v_' solid]);
rho_fluid = eval(['rho_' fluid]);
v_fluid = eval(['v_' fluid]);

% frequency band to integrate over
flow = find(freq>=.5E6,1);
fhigh = find(freq>=1.1E6,1);
freqband = freq(flow:fhigh);

% transmission coef
hmax = 800; % max thickness excursion in microns
h = (-hmax:0.5:hmax)*1E-6;
alpha = 2*pi*freqband*h/v_fluid;  % freq * thickness
Zr = rho_fluid*v_fluid/(rho_solid*v_solid);
rff = (Zr-1)/(Zr+1);
Trans = ((1-rff^2)*exp(-1i*alpha))./(1-rff^2*exp(-2*1i*alpha));

% objective function
Fun = zeros(nq,length(h));
% ref waveform
iiref = 3;
for ii = 1:nq
    Fun(ii,:) = sum(abs((U(flow:fhigh,ii)*ones(size(h))-(U(flow:fhigh,iiref)...
        *ones(size(h)).*Trans)).^2),1);
end

% locate min
[~, hmin] = min(Fun,[],2);
figure
disp('plotting fluid layer thikness vs time')
plot(AcqTime(1:end),h(hmin)*1E6)
xlabel('Acquisition time')
ylabel('Fluid layer thickness (\mum)')
datetick('x',15)
xlim([AcqTime(1) AcqTime(end)])
ylim([-0.5 1]*300)
title(['Source-receiver pair ' num2str(kk) '-' num2str(jj)])



% plots
figure
hold on
for ii = 226:228
    plot(h*1E6,Fun(ii,:),h(hmin(ii))*1E6,Fun(ii,hmin(ii)),'ok')
end
xlabel('Thickness (\mum)')
ylabel('Objective function')
legend({' num2str(ii)'})

%% nice figure for poster
% make figure
figure
set(gcf,'Position',[2286,1,1000,400])
yyaxis left
plot(AcqTime(1:end),h(hmin)*1E6)
datetick('x',15)
xlabel('Time')
ylabel('Fluid thickness (\mu m)')
xlim([PressureTime(1) PressureTime(end)])
yrange = [-0.2 1]*300;
ylim(yrange)
title(['Transmitted wave analysis for pair ' num2str(kk)])
hold on
plot(PressureTime([t1flow t1flow]),yrange,'k:')
%plot(PressureTime([t2flow t2flow]),yrange,'k:')
plot(PressureTime([tinit tinit]),yrange,'k:')
plot(PressureTime([tatm tatm]),yrange,'k:')

% axes('Position',[.34 .55 .25 .25])
% box on
% plot(AcqTime(tinit2-5:tinit2+7),h(hmin(tinit2-5:tinit2+7))*1E6)
% axis tight
% title('Breakdown inset')

%% add to other figure
hold on
yyaxis left
plot(AcqTime(1:end),h(hmin)*1E6,'-','Color',[0.3010 0.7450 0.9330])
linehandle = findobj(gca,'Type','line');

legend([linehandle(5) linehandle(1)],'pair 5','pair 6')
title('Transmitted wave analysis for pairs 5 and 6')

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
figure
disp('plotting fluid layer thickness vs time (2)')
plot(AcqTime(1:end),h(hmin)*1E6)
xlabel('Acquisition time')
ylabel('Fluid layer thickness (\mum)')
axis([datenum([AcqTime(1) AcqTime(end)]) [-2.5 0.5]*100])
legend(cellstr(num2str(SRpairs','%d')));

%% add events to plot
% define times // for 2018-05-28 data
% EventTimes = datetime({'14:22:00','14:36:00','14:45:00','14:50:00',...
%     '15:55:50','16:11:40','16:25:52','16:50:30'},'Format','HH:mm:ss');

% define times// for 2018-05-09 data
EventTimes = datetime({'10:00:00','10:10:00','10:14:00','10:30:00',...
   '11:30:00','11:55:00','12:22:00'},'Format','HH:mm:ss');

EventTimes.Year = 2000+datayear;
EventTimes.Month = datamonth;
EventTimes.Day = dataday;
% define vertical range
VertRange = [-1 1]*1E3;
hold on
axtmp = axis;
for ii = 1:length(EventTimes)
    plot([EventTimes(ii) EventTimes(ii)],VertRange,':k')
end
axis(axtmp)

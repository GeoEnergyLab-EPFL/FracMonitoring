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
endfile = 211411;

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
% updated aqcuisition VI from EPSLog and header info
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
nn = 5;

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


N = sqrt(sum(Dd.^2,1));

% check all on-axis pairs
figure
plot(T*1E6,Dd)

%% save selected pairs in new array
Pairs = {3,11,17,23};
DataS = zeros(nq,ns,length(Pairs));
for ii = 1:length(Pairs)
    DataS(:,:,ii) = squeeze(data3(:,:,Pairs{ii},Pairs{ii}));
end

%% pressure profile
fid = fopen([datafold 'voltage_' num2str(endfile) '.txt'],'r');
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

%% graphical display of amplitudes by source position for initial time step
% transducer positions
Pos = [zeros(1,8) 4:-1:1 -1:-1:-4; 4:-1:1 -1:-1:-4 zeros(1,8)]';

% max amplitude at arrival times
idx = find(T>=40E-6,1);
[M, Idx] = max(abs(data0(idx:end,:,:)));
M = squeeze(M);
Idx = squeeze(Idx);

% colors for max amplitude
Mnorm = M./max(max(M));
figure;
cmap = colormap(jet);
close gcf

d = 0.5; % marker diameter

if strcmp(debugmode,'on')
    % display transducers on map

    for kk = 1:nt
        figure
        for jj = 1:nr
            tmp_rect = rectangle('Position',[Pos(jj,1)-d/2,Pos(jj,2)-d/2,d,d],'Curvature',[1 1]);
            tmp_color = cmap(floor(length(cmap)*Mnorm(jj,kk)),:);
            set(tmp_rect,'FaceColor',tmp_color);
            if jj==kk
                set(tmp_rect, 'EdgeColor',[1 0 0]);
            end
        end
        axis equal
        axis([min(Pos(:,1))-d max(Pos(:,1))+d min(Pos(:,2))-d max(Pos(:,2))+d])
        title(['Source nb ' num2str(kk-1)])
        xlabel('x-position')
        ylabel('y-position')
        colormap(jet)
        colorbar
    end
end
%% test of distance normalization
% kk = 13; % reference source/receiver position
% Dist = ((Pos(:,1)-Pos(kk,1)).^2+((Pos(:,2)-Pos(kk,2)).^2)).^(1/2);
% Mdist = Mnorm(:,kk).*Dist/max(Mnorm(:,kk).*Dist);
% 
% % plot
% figure
% for jj = 1:nr
%     tmp_rect = rectangle('Position',[Pos(jj,1)-d/2,Pos(jj,2)-d/2,d,d],'Curvature',[1 1]);
%     if jj==kk
%         set(tmp_rect,'FaceColor','k');
%     else
%         tmp_color = cmap(ceil(length(cmap)*Mdist(jj)),:);
%         set(tmp_rect,'FaceColor',tmp_color);
%     end
% end
% axis equal
% axis([min(Pos(:,1))-d max(Pos(:,1))+d min(Pos(:,2))-d max(Pos(:,2))+d])
% title(['Source nb ' num2str(kk)])
% xlabel('x-position')
% ylabel('y-position')
% colormap(jet)
% colorbar
% 
% % figure to check relationship between distance and amplitude
% Mnorm2 = Mnorm(:,kk)+0.4.*Dist;
% 
% figure, hold on
% plot(Dist(1:4),Mnorm2(1:4),'ob')
% plot(Dist(5:8),Mnorm2(5:8),'oc')
% plot(Dist(9:12),Mnorm2(9:12),'or')
% plot(Dist(13:16),Mnorm2(13:16),'ok')
% xlabel('Distance (a.u.)')
% ylabel('Amplitude (a.u.)')

%% check source and receiver coupling
% generalize code for different source and receiver positionings
PosR = Pos; % receiver positions
PosT = Pos; % source positions
PosRg = repmat(PosR,nt,1);
PosTg = reshape(repmat(PosT,1,nr)',2,nr*nt)';
Dist = reshape(sqrt(sum((PosRg-PosTg).^2,2)),nr,nt);

% sum over all receivers for source strength
coef = 0.0;
Msrc = sum(Mnorm+coef*Dist,1)/max(sum(Mnorm+coef*Dist,1));
% plot
figure
for kk = 1:nt
    tmp_rect = rectangle('Position',[Pos(kk,1)-d/2,Pos(kk,2)-d/2,d,d],'Curvature',[1 1]);
    tmp_color = cmap(floor(length(cmap)*Msrc(kk)),:);
    set(tmp_rect,'FaceColor',tmp_color);
end
axis equal
axis([min(Pos(:,1))-d max(Pos(:,1))+d min(Pos(:,2))-d max(Pos(:,2))+d])
title('Source strength')
xlabel('x-position')
ylabel('y-position')
text(Pos(1,1)+0.5,Pos(1,2),'0')
text(Pos(8,1)+0.5,Pos(8,2),'7')
text(Pos(9,1)-0.25,Pos(9,2)+0.75,'8')
text(Pos(16,1)-0.25,Pos(16,2)+0.75,'15')
colormap(jet)
colorbar

% sum over all sources for receiver strength
coef = 0.0;
Mrec = sum(Mnorm+coef*Dist,2)/max(sum(Mnorm+coef*Dist,2));
% plot
figure
for ii = 1:nr
    tmp_rect = rectangle('Position',[Pos(ii,1)-d/2,Pos(ii,2)-d/2,d,d],'Curvature',[1 1]);
    tmp_color = cmap(floor(length(cmap)*Mrec(ii)),:);
    set(tmp_rect,'FaceColor',tmp_color);
end
axis equal
axis([min(Pos(:,1))-d max(Pos(:,1))+d min(Pos(:,2))-d max(Pos(:,2))+d])
title('Receiver strength')
xlabel('x-position')
ylabel('y-position')
text(Pos(1,1)+0.5,Pos(1,2),'0')
text(Pos(8,1)+0.5,Pos(8,2),'7')
text(Pos(9,1)-0.25,Pos(9,2)+0.75,'8')
text(Pos(16,1)-0.25,Pos(16,2)+0.75,'15')
colormap(jet)
colorbar

%% analyze waveform changes from fluid layer
% select source-receiver pair and corresponding data
jj = 1;
kk = 9; 
datatmp = squeeze(data3(:,:,jj,kk))';

% % quick figure
% figure
% plot(T*1E6,datatmp)
% axis([88 112 [-1 1]*5E-3])
% xlabel('Time (\mus)')
% ylabel('Amplitude (a.u.)')
% celllgd = cellstr(num2str(filenum','%4d'));
% legend(celllgd{:})

% lowpass filter
flowpass = 4E6; % cut at 5 MHz
[b, a] = butter(2,flowpass/Fn);
datatmp2 = filtfilt(b,a,datatmp);

% quick figure
figure
plot(T*1E6,datatmp2(:,1:3))
axis([40 80 [-1 1]*1])
%axis([88 112 [-1 1]*5E-1])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
celllgd = cellstr(AcqTime);
legend(celllgd{1:3})

% max as function of time
ampmax = max(abs(hilbert(datatmp2)));
figure, plot(AcqTime, ampmax,'o-')
axis([datenum([AcqTime(1) AcqTime(end)]) [0 8]*1E-1])
xlabel('Time')
ylabel('Max Amplitude')

% figure(3)
% plot(AcqTime, ampmax,'o-')

%% test issues with amplitude: missing signals before average?
correction = ones(size(ampmax));
correction([19 25 29 43 61 62 75 76 81]) = 50/49;
correction([2 66 99 114]) = 51/49;
correction(80) = (50/49);
ampmax2 = ampmax.*correction;
figure, plot(AcqTime, ampmax2,'o-')
axis([datenum([AcqTime(1) AcqTime(end)]) [0 8]*1E-1])
xlabel('Time')
ylabel('Max Amplitude')

%%
% % look at the last xx acquisitions
% xx = 4;
% 
% figure
% plot(T*1E6,datatmp2(:,end-xx+1:end))
% axis([88 112 [-1 1]*5E-3])
% xlabel('Time (\mus)')
% ylabel('Amplitude (a.u.)')
% celllgd = cellstr(num2str(filenum(end-xx+1:end)','%4d'));
% legend(celllgd{:})
% 
% % and at first ones
% xx = 8;
% figure
% plot(T*1E6,datatmp2(:,1:xx))
% axis([88 112 [-1 1]*5E-3])
% xlabel('Time (\mus)')
% ylabel('Amplitude (a.u.)')
% celllgd = cellstr(num2str(filenum(1:xx)','%4d'));
% legend(celllgd{:})

%% compare last ones to initial average
xx = 4;
% initial average and std dev
dataavg = mean(datatmp2(:,1:end-xx),2);
datastd = std(datatmp2(:,1:end-xx),0,2);

% figure with initial average and uncertainty
figure
colrs = get(gca,'colororder');
plot(T*1E6,dataavg,'LineWidth',2)
axis([88 112 [-1 1]*5E-3])
xlabel('Time (\mus)')
ylabel('Amplitude (a.u.)')
hold on
T2 = [T; flipud(T)];
areafill = [dataavg-2*datastd; flipud(dataavg+2*datastd)];
fill(T2*1E6, areafill,[0 0.447 0.741],'FaceAlpha',0.5,'LineWidth',0.5,'EdgeColor',colrs(1,:));

% add last four measurements
for yy = 1:xx
    plot(T*1E6,datatmp2(:,end-yy+1),'Color',colrs(yy+1,:))
end

%% freq domain for thickness estimation
% data with only inline source-receiver pairs
data4 = zeros(nq,ns,nr);
for ii = 1:nr
    data4(:,:,ii) = squeeze(data3(:,:,ii,ii));
end

% move to frequency domain
nfft = 2^nextpow2(ns);
U = fft(data4,nfft,2);
U = U(:,1:nfft/2+1,:)/nfft;
freq = Fn*linspace(0,1,nfft/2+1)';

% graphical check
kk = 5;
figure
plot(freq*1E-6,abs(U(:,:,kk)))%,freq*1E-6,abs(U(end,:,kk)))
xlabel('Frequency (MHz)')
ylabel('Amplitude (a.u.)')
axis([0 2 0 1E-2])

%%
% % look at the last xx acquisitions
% xx = 4;
% figure
% plot(freq*1E-6,abs(U(end-xx+1:end,:,kk)))
% axis([0 2 0 5E-5])
% xlabel('Frequency (MHz)')
% ylabel('Amplitude (a.u.)')
% celllgd = cellstr(num2str(filenum(end-xx+1:end)','%4d'));
% legend(celllgd{:})
% 
% % and at all excepted last four
% figure
% plot(freq*1E-6,abs(U(1:end-xx,:,kk)))
% axis([0 2 0 5E-5])
% xlabel('Frequency (MHz)')
% ylabel('Amplitude (a.u.)')
% celllgd = cellstr(num2str(filenum(1:end-xx)','%4d'));
% legend(celllgd{:})

%% thickness estimation
% experimental constants
rho_plexi = 1180;   % density PMMA
v_plexi = 2790;     % velocity PMMA
rho_glycerol = 1260;% density glycerol
v_glycerol = 1960;  % velocity glycerol

% frequency band to integrate over
flow = find(freq>=.6E6,1);
fhigh = find(freq>=.9E6,1);
freqband = freq(flow:fhigh);

% transmission coef
h = (-250:0.5:250)*1E-6;
alpha = 2*pi*freqband*h/v_glycerol;  % freq * thickness
Zr = rho_glycerol*v_glycerol/(rho_plexi*v_plexi);
rff = (Zr-1)/(Zr+1);
Trans = ((1-rff^2)*exp(-1i*alpha))./(1-rff^2*exp(-2*1i*alpha));

% objective function
kk = 5;
Fun = zeros(nq,length(h));
for ii = 1:nq
    Fun(ii,:) = sum(abs((U(ii,flow:fhigh,kk)'*ones(size(h))-(U(1,flow:fhigh,kk)'*ones(size(h)).*Trans)).^2),1);
end

% locate min
[~, hmin] = min(Fun,[],2);
figure, plot(AcqTime(1:end),h(hmin)*1E6,'o:')
xlabel('Time')
ylabel('Fluid layer thickness (\mum)')
axis([datenum([AcqTime(1) AcqTime(end)]) [-1 1]*10])

% plots
figure
plot(h*1E6,Fun(5,:))
xlabel('Thickness (\mum)')
ylabel('Objective function')
%axis([h(1) h(end) min(min(Fun)) max(max(Fun))])
% 
% figure
% plot(h*1E6,Fun(1:end-xx,:))
% xlabel('Thickness (\mum)')
% ylabel('Objective function')

%% Look at normal transmissions only
tt = 50;

datanorm = zeros(ns,nr);
for ii = 1:nr
    datanorm(:,ii) = squeeze(data3(tt,:,ii,ii));
end

% color image plot
figure
imagesc(0:nr-1,T*1E6,datanorm)
axis([0 nr-1 0 80])
caxis([-1 1]*0.5)
title('Normal transmission')
xlabel('Transducer number')
ylabel('Time (\mus)')
colormap(jet)

% waveform plot
figure
plot(T*1E6,datanorm)
axis([40 60 [-1 1]*1])
title('Normal transmission')
xlabel('Time (\mus)')
ylabel('Amplitude (V)')

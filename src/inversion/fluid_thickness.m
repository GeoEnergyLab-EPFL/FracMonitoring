function w = fluid_thickness(dataset,refsignal,T,...
    solidProperties,fluidProperties,timeWindow,freqWindow)

% function performing fluid thickness estimation on dataset with reference refsignals
% from Groenenboom & Fokkema (1998) https://doi.org/10.1190%2F1.1444306
%
% dataset indexed with the following format: sequences, time points, source-receiver pairs
% refsignals indexed with format: time points, source-receiver pairs
% T is the acquisition time vector
% Fn is the Nyquist frequency
% solidProperties is a structure containing the properties of the solid medium
% fluidProperties is a structure containing the properties of the fluid
% timeWindow is a 2-element vector with the limits for the window in time to consider
% freqWindow is a 2-element vector with the limits for the window in frequency to consider
%
% the function outputs w the computed width

np = length(T); % number of time points
Fs = 1/(T(2)-T(1)); % sampling frequency
Fn = Fs/2;  % Nyquist frequency
nq = size(dataset,1); % number of sequences
npairs = size(dataset,3); % number of source-receiver pairs to analyze

% create custom Hanning window
tmin = find(T>=timeWindow(1),1);
tmax = find(T>=timeWindow(end),1);
npHanning = 50;
winHanning = hanning(npHanning*2)';
Twindow = [zeros(1,tmin-npHanning) winHanning(1:npHanning) ones(1,tmax-tmin)...
    winHanning(npHanning+1:end) zeros(1,np-(tmax+npHanning))]';

% window signal with it
Twindow2D = repmat(Twindow,1,npairs);
Twindow3D = permute(repmat(Twindow,1,nq,npairs),[2,1,3]);
datasetWin = dataset.*Twindow3D;
refsignalWin = refsignal.*Twindow2D;

% move to freq domain
nfft = 2^nextpow2(np);
Udata = fft(datasetWin,nfft,2);
Udata = Udata(:,1:nfft/2+1,:)/nfft;
freq = Fn*linspace(0,1,nfft/2+1)';
% for ref signal
Uref = fft(refsignalWin,nfft);
Uref = Uref(1:nfft/2+1,:)/nfft;

% frequency band
flow = find(freq>=freqWindow(1),1);
fhigh = find(freq>=freqWindow(end),1);
freqband = freq(flow:fhigh);

% transmission coef
hmax = 60; % max thickness excursion in microns % 500 for MARB
% h = (-hmax:0.5:hmax)*1E-6; % vector of thicknesses
h = (-hmax:0.5:hmax)*1E-6; % vector of thicknesses
%h = (-0.01*hmax:0.5:hmax)*1E-6; % vector of thicknesses
alpha = 2*pi*freqband*h/fluidProperties.Vp;  % freq * thickness
Zr = fluidProperties.rho*fluidProperties.Vp/...
    (solidProperties.rho*solidProperties.Vp); % impedance contrast
rff = (Zr-1)/(Zr+1); % reflection coefficient
Trans = ((1-rff^2)*exp(-1i*alpha))./(1-rff^2*exp(-2*1i*alpha)); % transmission coefficient

% thin-layer assumption
% Trans = 2./(2+1i*alpha./Zr);

% objective function
Fun = zeros(nq,npairs,length(h));

% loop on acquisition sequences
for i_seq = 1:nq
    % now loop on source-receiver pairs
    for i_pair = 1:npairs
        % CAREFUL we take the NON-conjugate transpose here to get the right
        % dimensions
        Fun(i_seq,i_pair,:) = sum(abs((Udata(i_seq,flow:fhigh,i_pair).'*...
            ones(size(h))-(Uref(flow:fhigh,i_pair)*ones(size(h)).*Trans)).^2),1);
    end
end

% locate estimated thicknesses from min of objective function
[~, hmin] = min(Fun,[],3);
w = h(hmin);

% test plot

% for i_n=1:16
%     figure
%     plot(h,squeeze(Fun(40,i_n,:))) % for Local Seq.40
% end

end

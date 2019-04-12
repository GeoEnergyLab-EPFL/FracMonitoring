function [np, ns, nr, Fs, Fn, t0, dt, T] = read_acoustic_params(ActiveAcousticInfos)
% read active acoustic parameters from JSON header

% extract relevant parameters
np = ActiveAcousticInfos.NumberOfPoints;
ns = ActiveAcousticInfos.NumberOfSources;
nr = ActiveAcousticInfos.NumberOfReceivers;
Fs = ActiveAcousticInfos.SamplingFrequency_MHz_*1E6;
dt = 1/Fs;  % time step
t0 = 0;     % initial time
T = t0+dt*linspace(0,np-1,np)'; % time vector
Fn = 0.5*Fs;    % Nyquist frequency (Hz)

end
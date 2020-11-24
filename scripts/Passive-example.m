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

% 2020

datayear = 20;
%
% test on gabbro
datamonth = 11;
dataday = 04;
starttime = '141523';

% data folder name from experiment date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];

% extract header info from JSON file
fjson = [datapath datafold num2str(starttime) '.json'];
[jsonhdr,ActiveTransducers,myPlattens,myBlock,PassiveTransducers] = load_header(fjson);

%%

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

% read passive acoustic info from  JSON header
PassiveAcoustic = jsonhdr.PassiveAcousticInfos;



%%
xyz_passive = calc_global_coord(PassiveTransducers,myPlattens);

%% Plots passive sensor location

fig_b=plotblockwithplattens(myBlock,myPlattens);

fig_p=transducerplot3D(PassiveTransducers,myPlattens,fig_b);

%fig_a=transducerplot3D(myTransducers,myPlattens,fig_p);

 

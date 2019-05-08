
%% Function to measure sound velocity by reading two signals source-receiver at the begining


function [soundVelocity, soundVelocityPart, soundVelocityGI] = ...
    velocity_measurements(dataRef, dataset, Tmin, Tmax, dt, method, blockSizes)

% Ref signal:Swave = 162209, Pwave = 163814

% HOWTO call:
% EXAMPLE: [soundVelocity, soundVelocityPart, soundVelocityGI] = ... 
%     velocity_measurements('../data/19-01-11/163814', '../data/19-01-29/133723', 1, 1, 1, 42, 55)
% 
% (dataRefPath,dataSetPath,sequences,sources,receivers, Tmin, Tmax)
% dataRefPath -- path of the dataset, without an extantion
% dataSetPath -- path of the dataset, without an extantion
% sequences -- # of sequences in the dataSet
% sources -- # of source in the dataSet
% receivers -- # of receiver in the dataSet
% Tmin -- start of a timeframe for a direct signal, \mus
% Tmax -- end of a timeframe for a direct signal, \mus
% dt -- time step
% method -- which method do you want to start - Full data correlation, or selected Part or Graphical Input
% blockSizes -- size of Block, 3 directions

%Plot setting
set(0,'DefaultAxesColor',[1 1 1]) 
set(0,'DefaultLineLineWidth',3) %2
set(0,'DefaultAxesFontSize',24) %18
set(0,'DefaultTextFontSize',22) %14
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultTextFontWeight','bold')
set(0,'DefaultLineColor','k')

% ADDPATH
addpath('../src/data_processing/');
addpath('../src/data_load/');
addpath('../src/utilities/');
addpath('../src/properties/');
addpath('../src/inversion/');
addpath('../src/forward/');
addpath('../src/experiment_configuration/');


[~,nr] = size (dataset(1,:)); %number of receivers


%% Length between source and receiver
   
% length = 250; % normally from SampleBlock -- TBD
length = blockSizes(1);

% time step
% dt = timestep;

    
%% load data in dataSet
%%q dataSetFilePath = strcat(dataSetPath,'.bin')
%dataSet = cell(0,0) % need to be updated for cell
%for n = 1:16 % for now only 1 pair -- all others per each side TBD
%%q    dataSet{1} = load_data(dataSetFilePath,sequences,sources,receivers); %data path !!! N !!!
%end
    dataSet = dataset(:,1); %data path !!! N !!!


%% reference signal read (without a sample)
%%q dataRefFilePath = strcat(dataRefPath,'.bin')
%%q dataRef = load_data(dataRefFilePath,1,1,1);

% dataRef = dataRef;


% figure % plot dataRef
% disp('reference signal read, source')
% plot(dataTime*1E6,dataRef,dataTime*1E6,abs(hilbert(dataRef)))


% figure % plot dataSet
% disp('plot for dataSet, receiver')
% plot(dataTime*1E6,dataSet)


%% Find t_propagation and Velocity with xcorr dataSet and dataRef
if (find(contains(method, 'Full')) > 0)

    [acor, lag] = xcorr(dataSet, dataRef);

    % % Figure check xcorr to lag
    % figure
    % disp('Xcor plot')
    % plot (lag*dt,acor)

    %t_propagation
    [~,i] = max(acor); % find max of correlation, i = position
    t_propagation = lag(i)*dt; % lagTime = time propagation

    % Velocity
    soundVelocity = length / t_propagation % m/s

else
    soundVelocity = 'Not selected';
end


%% PART: Find t_propagation and Velocity with xcorr dataSetPart and dataRefPart
if (find(contains(method, 'Part')) > 0)

    % set Tmin and Tmax
    Tmin = Tmin*1E-6/dt; % SET Tmin for dataSet <--
    Tmax = Tmax*1E-6/dt; % SET Tmax for dataSet <--

    % take part of dataSet and dataRef
    dataRefPart = dataRef(1:1000-1,:);
    %%q dataSetPart = dataSet{1,1}; 
    dataSetPart = dataSet(round(Tmin):round(Tmax),:);


    % PLOT dataSet + dataSetPart on 1 figure
    % figure % plot dataSet
    % plot(dataTime*1E6,dataSet{1})
    % hold on
    % plot(dataTime(:,round(Tmin):round(Tmax))*1E6,dataSetPart,'.r')
    % legend('dataSet','dataSetPart')
    % title('dataSet + dataSetPart')


    % figure % plot dataRefPart
    % disp('Xcor part, dataRefPart')
    % plot(dataTime(:,1:1000-1),dataRefPart)

    % figure % plot dataSetPart
    % plot(dataTime(:,Tmin:Tmax),dataSetPart)


    % xcorr with dataRefPart and dataSetPart
    [acorPart, lagPart] = xcorr(dataSetPart, dataRefPart);

    %Figure check xcorr to lag
    % figure 
    % disp('Xcor plot, PART')
    % plot (lagPart*dt,acorPart)

    %t_propagation
    [~,j] = max(acorPart); % find max of correlation, i = position
    %lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
    t_propagationPart = (Tmin+lagPart(j))*dt; % lagTime = timepropagation

    % VelocityPart
    soundVelocityPart = length / t_propagationPart % m/s

else
    soundVelocityPart = 'Not selected';
end


%% Check with GI input
if (find(contains(method, 'GI')) > 0)
    
    
    %take time step (dt) for data -- normally header -- TBD
    % dataSetJsonPath = strcat(dataSetPath,'.json')
    % [~, ~, jsonhdr] = load_header(dataSetJsonPath); % read from header
    % Fs = jsonhdr.ActiveAcousticInfos.SamplingFrequency_MHz_*1E6; %freq of measurements 
    % dt = 1/Fs;  % time step
    % dt = timestep;

    %create a timelist for easy plotting
    dataTime(:,1) = 0;
    for i = 1:7999
        dataTime(:,i+1) = dataTime(:,i)+dt;
    end
    
    figure(gcf)
    disp('GI input')
    hold on
    plot(dataTime*1E6,dataSet*10^4) %dataSet
    plot(dataTime*1E6,dataRef*10^3) %dataRef
    plot(dataTime(:,round(Tmin):round(Tmax))*1E6,dataSetPart*10^4,'.g') %dataPart
    xlim([0,100]) %100/160 by X <---
    xlabel('Time, \mus')
    ylabel('Amplitude, mV')
    legend('Receiver signal x10 Scale!','Source signal', 'dataSetPart for xcorr')
    title('Amlitude of Source and Receiver')


    [x,y] = ginput(2)
    t_propagationGI = abs(x(2) - x(1))
    soundVelocityGI = length * 10^6 / t_propagationGI % m / timeScale
    % soundVelocityGI = 0;
else
    soundVelocityGI = 'Not selected';
end


end


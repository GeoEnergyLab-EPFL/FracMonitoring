
%% Function to measure sound velocity by reading two signals source-receiver at the begining


function [soundVelocity, soundVelocityPart, soundVelocityGI, soundVelocityAverage] = ...
    velocity_measurements(dataRef, dataSet, Tmin, Tmax, method, dt, myBlock)

% Ref signal:Swave = 162209, Pwave = 163814

% HOWTO call:
% EXAMPLE: [Vp, VpPart, VpGI, VpMedian] = velocity_measurements(dataRef, dataset, 35, 45, ["Full", "Part", "GI"], ... % dataRef, dataset, Tmin, Tmax, ["Full", "Part", "GI"], 
%    dt, myBlock.sizes, jsonhdr.ActiveAcousticInfos.ArrayPiezoPosition); % data taken from results above 
% 
% dataRef -- path of the dataset for Ref signal (source)
% dataset -- path of the dataset for receiver
% Tmin -- 35, start of a timeframe for a direct signal, \mus
% Tmax -- 45, end of a timeframe for a direct signal, \mus
% method -- ["Full", "Part", "GI"], which method do you want to start - Full data correlation, or selected Part or Graphical Input
% dt -- time step
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
addpath('../FracMonitoring/src/data_processing/');
addpath('../FracMonitoring/src/data_load/');
addpath('../FracMonitoring/src/utilities/');
addpath('../FracMonitoring/src/properties/');
addpath('../FracMonitoring/src/inversion/');
addpath('../FracMonitoring/src/forward/');
addpath('../FracMonitoring/src/experiment_configuration/');


        dataTime(:,1) = 0;
        for i = 1:7999
            dataTime(:,i+1) = dataTime(:,i)+dt;
        end

 figure % plot dataRef
 disp('reference signal read, source')
 plot(dataTime*1E6,dataRef)

[~,nr] = size (dataSet(1,:)); %number of receivers


%% Length between source and receiver
   
% length = 250; % normally from SampleBlock -- TBD
length = myBlock.sizes(1);


%% Find t_propagation and Velocity with xcorr dataSet and dataRef
if (find(contains(method, 'Full')) > 0)
    for ii = 1:nr
        [acor, lag] = xcorr(dataSet(:,ii), dataRef);

        % % Figure check xcorr to lag
        % figure
        % disp('Xcor plot')
        % plot (lag*dt,acor)

        %t_propagation
        [~,i] = max(acor); % find max of correlation, i = position
        t_propagation = lag(i)*dt; % lagTime = time propagation

        % Velocity
        soundVelocity(ii) = length / t_propagation % m/s
    end
    
else
        soundVelocity = 'Not selected';

end

Vp_median = median(soundVelocity);


%% PART: Find t_propagation and Velocity with xcorr dataSetPart and dataRefPart
if (find(contains(method, 'Part')) > 0)

        % set Tmin and Tmax
        Tmin = Tmin*1E-6/dt; % SET Tmin for dataSet <--
        Tmax = Tmax*1E-6/dt; % SET Tmax for dataSet <--

        % take part of dataSet and dataRef
        dataRefPart = dataRef(1:1000-1,:);
        %%q dataSetPart = dataSet{1,1}; 
        for ii = 1:nr
            dataSetPart(:,ii) = dataSet(round(Tmin):round(Tmax),ii);
        end


    
    for ii = 1:nr

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
        [acorPart, lagPart] = xcorr(dataSetPart(ii,:), dataRefPart);

         %figure check xcorr to lag
         %figure 
         %disp('Xcor plot, PART')
         %plot (lagPart*dt,acorPart)

        %t_propagation
        [~,j] = max(acorPart); % find max of correlation, i = position
        %lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
        t_propagationPart = (Tmin+lagPart(j))*dt; % lagTime = timepropagation

        % VelocityPart
        soundVelocityPart(ii) = length / t_propagationPart % m/s
    end
else
    soundVelocityPart = 'Not selected';

end

VpPart_median = median(soundVelocityPart);

%% Check with GI input
if (find(contains(method, 'GI')) > 0)
    
    
            %create a timelist for easy plotting
        dataTime(:,1) = 0;
        for i = 1:7999
            dataTime(:,i+1) = dataTime(:,i)+dt;
        end
        
    for ii = 1:nr 
        
        figure % (gcf)
        disp('GI input')
        hold on
        plot(dataTime*1E6,dataSet(:,ii)*10^3) %dataSet
        plot(dataTime*1E6,dataRef*10^3) %dataRef
        plot(dataTime(:,round(Tmin):round(Tmax))*1E6,dataSetPart(:,ii)*10^3,'.g') %dataPart
        xlim([0,100]) %100/160 by X <---
        xlabel('Time, \mus')
        ylabel('Amplitude, mV')
        legend('Receiver signal','Reference signal', 'Selected data')
        set(gca,'box','on')
        % title('Amlitude of Source and Receiver')



        [x,y] = ginput(2)
        t_propagationGI = abs(x(2) - x(1))
        soundVelocityGI(ii) = length * 10^6 / t_propagationGI % m / timeScale
        % soundVelocityGI = 0;

    end
else
    soundVelocityGI = 'Not selected';
end

VpGI_median = median(soundVelocityGI);

soundVelocityAverage = [Vp_median; VpPart_median; VpGI_median] % 101 means no value

end




% DO NOT LOOK BELOW


%% load data in dataSet
%%q dataSetFilePath = strcat(dataSetPath,'.bin')
%dataSet = cell(0,0) % need to be updated for cell
%for n = 1:16 % for now only 1 pair -- all others per each side TBD
%%q    dataSet{1} = load_data(dataSetFilePath,sequences,sources,receivers); %data path !!! N !!!
%end
%    dataSet = dataset(:,1); %data path !!! N !!!


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

close all
clearvars
home

import SourceReceiver

Pair = [];
for i = 1:32:1
    Pair = Pair + SourceReceiver(get_mesh_coord(S(i,1),S(i,2),S(i,3)),get_mesh_coord(R(i,1),R(i,2),R(i,3)));
end

%% extract coordinates on ENACDrives

if isunix
    [~, uid] = system('id -u');
    [~, username] = system('whoami');
    coordpath = ['/run/user/' uid(1:end-1) '/gvfs/smb-share:domain=INTRANET,server=enac1files.epfl.ch,'...
     'share=gel,user=' username(1:end-1) '/research/Experiments/HF-Experiments/GEL-data/'];
 elseif ispc
     datapath = 'Z:/research/Experiments/HF-Experiments/GEL-data/';
end

fid = fopen('source_receivers_coordinates.txt', 'rt');
nlines = 0;
while (fgets(fid) ~= -1 )
    nlines = nlines+1;
end    
fclose(fid);

S1 = textscan('source_receivers_coordinates', '%s','delimiter', '\n'); % http://matlab.izmiran.ru/help/techdoc/ref/textscan.html
lm1 = 2;
rowindx = 1;

S2 = textscan('source_receivers_coordinates', '%s','delimiter', '\n');
lm1 = 3;
rowindx = 1;




%% defining the constants

d = 0.25; % size of the sample, m
n_x = 101; % mesh of x
n_y = 101; % mesh of y
n_z = 101; % mesh of z
nt = 32;  % number of sources
nr = 32;  % number of receivers
Fs = 5E7;   % sampling frequency (hardware defined)
dt = 1/Fs;  % time step
t0 = 0;     % initial time
T = t0+dt*linspace(0,ns-1,ns)'; % time vector
Fn = 0.5*Fs; % Nyquist frequency (Hz) (the minimum rate at which a signal can be sampled without introducing errors, which is twice the highest frequency present in the signal)

S1 = [10 d/2 10]; % coordintes of sources
S2 = [];
S3 = [];
S4 = [];

S = [S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24 S25 S26 S27 S28 S29 S30 S31 S32];
R = [R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28 R29 R30 R31 R32];


%% creation of the mesh

function [L_x,L_y,L_z] = create_mesh(n_x, n_y, n_z)
    for i = 1:n_x:1 %i = 0; i<n_x; i++  %from 1 to 100
        L_x(i) = (-d/2 + (i-1)*2*d/(n_x-1));
    end
    for i = 1:n_y:1  
        L_y(i) = (-d/2 + (i-1)*2*d/(n_y-1));
    end
    for i = 1:n_z:1  
        L_z(i) = (-d/2 + (i-1)*2*d/(n_z-1));
    end
end

%% assignment of coordinates

function [x,y,z] = get_mesh_coord(Measurementx,Measurementy,Measurementz)
        for i = 1:n_x:1
            if (i<n_x && Measurementx > L_x(i) && Measurementx < L_x(i+1))
                x = L_x(i);
            
            else 
                disp("error on x coordinate")
            end
        end
        for i = 1:n_y:1
            if (i<n_y && Measurementy > L_y(i) && Measurementy < L_y(i+1))
                y = L_y(i);
            
            else 
                disp("error on y coordinate")
            end
        end
        for i = 1:n_z:1
            if (i<n_z && Measurementz > L_z(i) && Measurementz < L_z(i+1))
                z = L_z(i);
            
            else 
                disp("error on z coordinate")
            end
        end
end


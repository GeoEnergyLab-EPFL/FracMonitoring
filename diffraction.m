close all
clearvars
home

import SourceReceiver

Pair = [];
for i = 1:32:1
    Pair = Pair + SourceReceiver(get_mesh_coord(S(i,1),S(i,2),S(i,3)),get_mesh_coord(R(i,1),R(i,2),R(i,3)));
end

Sources = [];
for i = 1:length(Pair):1
    Sources = Sources + Pair(i).get_source();
end

Receivers = [];
for i = 1:length(Pair):1
    Receivers = Receivers + Pair(i).get_receiver();
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
c = 5120; % speed of sound in sample (m/s) (here for iron)
S1 = [10 d/2 10]; % measured coordinates of sources
S2 = [];
S3 = [];
S4 = [];

S = [S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24 S25 S26 S27 S28 S29 S30 S31 S32];
R = [R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28 R29 R30 R31 R32];
mesh = create_mesh(n_x, n_y, n_z);

%% creation of the mesh

function [L_x,L_y,L_z] = create_mesh(n_x, n_y, n_z)
    for i = 1:n_x:1 
        L_x(i) = (-d/2 + (i-1)*2*d/(n_x-1));
    end
    for i = 1:n_y:1  
        L_y(i) = (-d/2 + (i-1)*2*d/(n_y-1));
    end
    for i = 1:n_z:1  
        L_z(i) = (-d/2 + (i-1)*2*d/(n_z-1));
    end
end

%% assignment of coordinates in meters

function [x,y,z] = get_mesh_coord(Measurementx,Measurementy,Measurementz)
        for i = 1:n_x:1
            if (i<n_x && Measurementx > mesh(1,i) && Measurementx < mesh(1,i+1))
                x = mesh(1,i); 
              
            end
            
        end
        for i = 1:n_y:1
            if (i<n_y && Measurementy > mesh(2,i) && Measurementy < mesh(2,i+1))
                y = mesh(2,i);
           
            end
        end
        for i = 1:n_z:1
            if (i<n_z && Measurementz > mesh(3,i) && Measurementz < mesh(3,i))
                z = mesh(3,i);
            end
        end
end

%% Generation of fracture given the major axis and the minor axis in meters
function Fracture = Get_Fracture(LengthMajorAxis,LengthMinorAxis)
    x = 0;
    while x < LengthMajorAxis
        y2 = LengthMinorAxis^2*(1 - x^2/LengthMajorAxis^2);
        i=0;
        if y2>=0
            y = sqrt(y2);
            for i = 1:n_y:1
                if (i<n_y && Measurementy > mesh(2,i) && Measurementy < mesh(2,i+1))
                    y = mesh(2,i);
                    
                end
            end
            Fracture= Fracture + [x,y];
            i = i + 1;
            x = mesh(1,i);
        else
            break
        end
    end
end

%% Function calculating the distance bewteen two points
function dist = Distance(A,B) 
    dist = sqrt((B(1)-A(1))^2 + (B(2)-A(2))^2);
end

%% Given a point  S, returns the closest point in a list of coordinates "Fracture"
function [Point,dist] = ShortestDistance(Source,Fracture)
    dist = Distance(Source,Fracture(1));
    Point = Fracture(1);
    for i = 2:length(Fracture):1
        dist2 = Distance(Source,Fracture(i));
        if dist2<dist
            dist = dist2;
            Point = Fracture(i);
        end
    end
end 

%% Returns the closest receiver to a given source after diffraction and the total distance (source-fracture + fracture-receiver)
function [ClosestR,DiffractedDist] = ClosestReceiver(Source,Fracture)
    res1 = ShortestDistance(Source,Fracture);
    res2 = ShortestDistance(res1(2),Receivers);
    DiffractedDist = res2(1)+res(1);
    ClosestR = res2(2);
end

%% Returns the timetable of all DeltaT for each source.
function Timetable = DeltaT(Fracture)
    for i = 1:length(Sources):1
        ReceiverResult = ClosestReceiver(Sources(i),Fracture);
        DirectTime = Distance(Sources(i),ReceiverResult(1))/c;
        DiffractedTime = ReceiverResult(2)/c;
        Timetable = Timetable + abs(DiffractedTime-DirectTime);
    end
end

function [timesolution, a, b] = Timesolutions()
    for i = 1:n_x:1
        for j = 1:n_y:1
            a = mesh(1,i);
            b = mesh(2,j);
            Fracture = Get_Fracture(a,b);
            timesolution = timesolution + DeltaT(Fracture);
        end
    end
end


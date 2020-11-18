function [json_header,activeTransd,myPlattens,myBlock,varargout] = load_header(filename)
%  
% purpose - load json header and create block objects, plattens &
% transducers objects
%
% last modification : Novmber 18 2020 - output also the passive
% transducers.

% general header info
json_header = jsondecode(fileread(filename));

% block info converted to meters
length_E = json_header.TestInfos.BlockDimensions_mm_.E_W_mm_*1E-3;
length_N = json_header.TestInfos.BlockDimensions_mm_.N_S_mm_*1E-3;
length_T = json_header.TestInfos.BlockDimensions_mm_.T_B_mm_*1E-3;
for direc = ['E','N','T']
    if isempty(eval(['length_' direc]))
        eval(['length_' direc ' = 250;'])
    end
end
myBlock = SampleBlock([length_E, length_N, length_T]);

%%%% PLATTEN
% empty platten array
n_platten = length(json_header.ActiveAcousticInfos.ArrayPlattenPosition);
myPlattens = Platten.empty(n_platten,0);
for i_platten = 1:n_platten
    % corrected bug - missing indices -> as a result  the first platten was
    % duplicated 6 times;(
    % 
    tmpID = json_header.ActiveAcousticInfos.ArrayPlattenPosition(i_platten).Platten;
    tmpFace = json_header.ActiveAcousticInfos.ArrayPlattenPosition(i_platten).Position;
    tmpPlatten = Platten(tmpID,tmpFace,myBlock);
    myPlattens(i_platten) = tmpPlatten;
end

%%%% ACTIVE TRANSDUCERS
% active acoustics transducer info
ns = json_header.ActiveAcousticInfos.NumberOfSources;
nr = json_header.ActiveAcousticInfos.NumberOfReceivers;
n_transd = ns+nr;

% create empty arrays for info
serial = zeros(n_transd,1);
type = blanks(n_transd)';
channel = zeros(n_transd,1);
platten = blanks(n_transd)';
local_id = zeros(n_transd,1);
orientation = zeros(n_transd,1);

% get all the active transducer info one by one
for i_transd = 1:n_transd
    transd_info = json_header.ActiveAcousticInfos.ArrayPiezoPosition(i_transd);
    serial(i_transd) = str2double(transd_info.PiezoSN);
    type(i_transd) = transd_info.Type;
    channel(i_transd) = str2double(transd_info.Channel);
    platten(i_transd) = transd_info.Platten;
    local_id(i_transd) = str2double(transd_info.Position);
    orientation(i_transd) = str2num(transd_info.Orientation); % keep str2num to properly read 'pi'
end
% create active transducers object
activeTransd = Transducers(serial,type,channel,platten,local_id,orientation);


%%%%%% PASSIVE TRANSDUCERS
%
n_transd_p = length(json_header.PassiveAcousticInfos.ArraySensorPositions);

 if  n_transd_p>0
% create empty arrays for info
serial = zeros(n_transd_p,1);
type = blanks(n_transd_p)';
channel = zeros(n_transd_p,1);
platen = blanks(n_transd_p)';
local_id = zeros(n_transd_p,1);
orientation = zeros(n_transd_p,1);

for i = 1:n_transd_p
    transd_info = json_header.PassiveAcousticInfos.ArraySensorPositions(i);
    serial(i) = str2double(transd_info.PiezoSN);
    type(i) = 'R';
    channel(i) = str2double(transd_info.Channel);
    platen(i) = transd_info.Platen;
    local_id(i) = str2double(transd_info.Position);
    orientation(i) = 0; % keep str2num to properly read 'pi'
end
varargout = cell(1,1);
% create active transducers object
if length(varargout)==1
    varargout{1} = Transducers(serial,type,channel,platen,local_id,orientation);
 end
 end
 
end

function [json_header,myTransd,myPlattens,myBlock] = load_header(filename)
% TODO: write code
% purpose - load json header and create block objects, plattens &
% transducers objects

% general header info
json_header = jsondecode(fileread(filename));

% block info
length_E = json_header.TestInfos.BlockDimensions_mm_.E_W_mm_;
length_N = json_header.TestInfos.BlockDimensions_mm_.N_S_mm_;
length_T = json_header.TestInfos.BlockDimensions_mm_.T_B_mm_;
myBlock = SampleBlock([length_E, length_N, length_T]);

%%%% PLATTEN
% empty platten array
n_platten = length(json_header.ActiveAcousticInfos.ArrayPlattenPosition);
myPlattens = Platten.empty(n_platten,0);
for i_platten = 1:n_platten
    tmpID = json_header.ActiveAcousticInfos.ArrayPlattenPosition.Platten;
    tmpFace = json_header.ActiveAcousticInfos.ArrayPlattenPosition.Position;
    tmpPlatten = Platten(tmpID,tmpFace,myBlock);
    myPlattens(i_platten) = tmpPlatten;
end

%%%% TRANSDUCERS
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

% get all transducer info one by one
for i_transd = 1:n_transd
    transd_info = json_header.ActiveAcousticInfos.ArrayPiezoPosition(i_transd);
    serial(i_transd) = str2double(transd_info.PiezoSN);
    type(i_transd) = transd_info.Type;
    channel(i_transd) = str2double(transd_info.Channel);
    platten(i_transd) = transd_info.Platten;
    local_id(i_transd) = str2double(transd_info.Position);
    orientation(i_transd) = str2num(transd_info.Orientation); % keep str2num to properly read 'pi'
end
% create transducers object
myTransd = Transducers(serial,type,channel,platten,local_id,orientation);

end

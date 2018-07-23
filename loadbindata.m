function [ dataInit2,dataInit3 ] = loadbindata(ns,nr,nt,nq,q0,qq,datapath2,FolderInfo,startindex,a,b)
%LOADBINDATA loads selected binary data files for QQ acquisition snapshots
%   starting at Q0, with parameters NS, NR, NT and NQ, and DC filter
%   coefficients A and B
%
% Thomas Blum, Geo-Energy Lab, EPFL, July 2018

% Error message if trying to load too many data files
if q0+qq-1>nq
    error('loadbindata:dataindex','Trying to load data sequence too far in time')
end

% create temporary arrays
datatmp = zeros(qq,ns,nr*nt);
datafilt = zeros(size(datatmp));

% read data and load it to temp arrays
for ii = q0:q0+qq-1
    fid = fopen([datapath2 FolderInfo(startindex+ii-1).name],'r');
    datatmp(ii-q0+1,:,:) = fread(fid,[ns,nt*nr],'double');
    fclose(fid);
    datafilt(ii-q0+1,:,:) = filtfilt(b,a,squeeze(datatmp(ii-q0+1,:,:))); % DC filter with a and b coefficients
end

% reshape data for outputs
dataInit2 = reshape(datatmp,qq,ns,nr,nt);
dataInit3 = reshape(datafilt,qq,ns,nr,nt);
% qq number of acquisition snapshots
% ns number of samples
% nr number of receivers
% nt number of sources

end
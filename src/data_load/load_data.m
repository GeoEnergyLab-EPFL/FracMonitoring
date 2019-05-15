function  [dataout, varargout] = load_data(filename,sequences,sources,receivers,varargin)
% function to load a set of sequences, sources, receivers data

% in varargin min-max time to reduce the number of data pts per traces (NOT
% YET CODED)
% return a 4D array
%  d-1 sequences
%  d-2 time
%  d-3 source
%  d-4 receivers

% check if header file exists
fjson =  strrep(filename,'.bin','.json');
if ~isfile(fjson)
        disp('no corresponding header file found')
    return
end
% get active acoustic info from header
jsonhdr = load_header(fjson);
np = jsonhdr.ActiveAcousticInfos.NumberOfPoints;
ns = jsonhdr.ActiveAcousticInfos.NumberOfSources;
nr = jsonhdr.ActiveAcousticInfos.NumberOfReceivers;

% size in bytes of one single acquisition sequence
seqsize = ns*nr*np*8;

% open bin data file and read header size (recorded in first double in bin)
fid = fopen(filename,'r');
hdrsize = (fread(fid,1,'double')+1)*8;

% check bin file size and get number of acquisition sequences
binstats = dir(filename);
if ~isequal(mod(binstats.bytes,seqsize),hdrsize)
    disp('incomplete acquisition sequence')
end

nq = floor(binstats.bytes/seqsize); % number of total acquisition sequences
qq = length(sequences); % number of acq sequences to load

% check that all sequences exist
if sum(sequences>nq)
    disp('trying to read non-existant sequences')
    return
end

% set pointer to beginning of binary data for q0 sequence
fseek(fid,hdrsize+(sequences(1)-1)*seqsize,'bof');

% create empty array
datatmp2 = zeros(qq,np,length(receivers),length(sources));
% display progress
textprogressbar('loading sequences:    ');
for ii = 1:length(sequences)
    % set pointer to beginning of sequence
    fseek(fid,hdrsize+(sequences(ii)-1)*seqsize,'bof');
    datatmp1 = fread(fid,[np,nr*ns],'double');
    datatmp1 = reshape(datatmp1,np,nr,ns);
    datatmp2(ii,:,:,:) = datatmp1(:,receivers,sources);
    pgrs = floor(ii/length(sequences)*100);
    textprogressbar(pgrs);
end
textprogressbar(' done');

fclose(fid);
% reshape data array to 4D with dimensions in intuitive order
dataout = permute(datatmp2,[1 2 4 3]);
varargout = {nq};
end

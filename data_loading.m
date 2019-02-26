function  dataout = data_loading(filename,sequences,sources,receivers,np,nr,ns,varargin)
% function to load a set of sequences, sources, receivers data

% in varargin min-max time to reduce the number of data pts per traces (NOT
% YET CODED)
% return a 4D array
%  d-1 sequences
%  d-2 time
%  d-3 source
%  d-4 receivers

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
datatmp = zeros(qq,np,nr*ns);
for ii = 1:length(sequences)
    % set pointer to beginning of  sequence
    fseek(fid,hdrsize+(sequences(ii)-1)*seqsize,'bof');
    datatmp(ii,:,:) = fread(fid,[np,ns*nr],'double');
end
fclose(fid);
% reshape data array to 4D with dimensions in intuitive order
datatmp2 = permute(reshape(datatmp,qq,np,nr,ns),[1 2 4 3]);
dataout = datatmp2(:,:,sources,receivers);

end

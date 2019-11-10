function  [Pair_info, Pair_acqT] = load_diffraction(fpath, sidemarker, wave_type, seqnb)
% Dong Liu -- 10/10/2019
% function to load the arrival time for a given sequence
%
% Input:
% fpath: path of the diffraction data
% sidemarker: the side you are interested in 'N-E-S-W'
% wave_type: the wave type you are interested in
% seqnb: the selected sequence number
%
% Output:
% Pair_info: [S channel number, R channel number and arrival time]
% Pair_acqT: acquisition time, format_date time
% count the pairs in all these four directions

[n_pair,full_path, SRmap] = getPairInfo(fpath, sidemarker, wave_type);

% get the arrival time and the source-receiver pair where the arrival time
% of the selected sequence is picked
Pair_info = zeros(sum(n_pair),3);
Pair_acqT = [];
i_ct = 0;
for i = 1:sum(n_pair)
    [arrival_info] = importdata(full_path{i},' ');
    % this will import the data into textdata(date and time) and data (sequence number, arrival time and source number, receiver number)
    picked_seqs = arrival_info.data(1:end,1);
    picked_s = arrival_info.data(1,3);
    picked_r = arrival_info.data(1,4);
    % picking for gabbro, the local sequence = the global sequence
    % one needs to verify again with the recording time, to be coded 
    [~, picked_idx] = ismember(seqnb,picked_seqs);
    if picked_idx > 0 && (picked_s == SRmap(i,1)) && (picked_r == SRmap(i,2))
        i_ct = i_ct+1;
        Pair_info(i_ct,:) = [picked_s picked_r arrival_info.data(picked_idx,2)];% the arrival time is in \mus
        if isempty(Pair_acqT)
            picked_seq = string(arrival_info.textdata(picked_idx,1:end));
            %disp(picked_seq);
            picked_seq = strcat(picked_seq(1:end,1)," ", picked_seq(1:end,2));
            Pair_acqT = datetime(picked_seq,'InputFormat','dd-MMM-yy HH:mm:ss');
        end
    end
end
Pair_info(i_ct+1:sum(n_pair),:) = []; % delete the rows where there is no information

% one sequence corresponding to one acquisition time

end


function [nPair,full_path, SRmap]=getPairInfo(fpath, sidemarker, wave_type)
    nPair = zeros(length(sidemarker),1);
    for i = 1:length(sidemarker)
        subfold = [fpath sidemarker(i,:) '/*' wave_type '.txt'];
        localnames = dir(subfold);
        nPair(i,1) = length(localnames);
    end
    total_pair=sum(nPair);
    % get the source and receiver channel
    full_path = cell(total_pair,1);
    SRmap = zeros(total_pair,2);
    pair_ct = 0;
    for i = 1:length(sidemarker)
        subfold = [fpath sidemarker(i,:) '/*' wave_type '.txt'];
        localnames = dir(subfold);
        for j = 1:length(localnames)
            pair_ct = pair_ct+1;
            [full_path{pair_ct},SRmap(pair_ct,1:2)] = read_names(localnames(j));
        end
    end
end


function[full_path,SR_map]=read_names(NameObj)
% Dong Liu -- 04/10/2019
% return the SR map knowing the directory object of the picked data and the
% full directory for one pair
full_path=[NameObj.folder '/' NameObj.name];
SR_map=readname(NameObj.name);
end


function[SR]=readname(name)
% return the S channel and R channel from the name
% the file ends with wave_type('PP' for example)+'.txt'
% After removing the last 6 characters:
% the S channel is from 2nd to the position before R
% the R channel is from the R position+1 to the end
real_name=name(1:end-6);
i_R=strfind(real_name,'R');
i_R=i_R(1);
SR=[str2num(real_name(2:i_R-1)) str2num(real_name(i_R+1:end))];

end

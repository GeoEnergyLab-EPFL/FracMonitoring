function [datapath] = pathbyarchitecture()
%PATHBYARCHITECTURE builds path depending on computer architecture
%   PATHBYARCHITECTURE returns the ENACdrives GEL-data path corresponding
%   to the computer architecture on which the function is running
%
% Thomas Blum, Geo-Energy Lab, EPFL, June 2018

if ismac % path on macos
    [~, username] = system('whoami');
    datapath = ['/Users/' username(1:end-1) '/ENACdrives/gel_on_enac1files/'...
        'research/Experiments/HF-Experiments/GEL-data/'];
elseif isunix % path on unix/linus os
    [~, username] = system('whoami');
    [~, hostname] = system('hostname');
    if strcmp(hostname(1:end-1),'X1Carbon') % special case for mount.cifs on X1Carbon
        datapath = ['/home/' username(1:end-1) '/ENACdrives/gel_on_enac1files/research/'...
            'Experiments/HF-Experiments/GEL-data/'];
    else % path for more typical GVFS mount
        [~, uid] = system('id -u');
        datapath = ['/run/user/' uid(1:end-1) '/gvfs/smb-share:domain=INTRANET,server=enac1files.epfl.ch,'...
            'share=gel,user=' username(1:end-1) '/research/Experiments/HF-Experiments/GEL-data/'];
    end
elseif ispc % path on windows
    datapath = 'Y:/research/Experiments/HF-Experiments/GEL-data/';
else
    disp('Platform not supported')
end
end


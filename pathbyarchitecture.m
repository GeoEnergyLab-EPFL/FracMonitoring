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
    [~, uid] = system('id -u');
    [~, username] = system('whoami');
    datapath = ['/run/user/' uid(1:end-1) '/gvfs/smb-share:domain=INTRANET,server=enac1files.epfl.ch,'...
        'share=gel,user=' username(1:end-1) '/research/Experiments/HF-Experiments/GEL-data/'];
elseif ispc % path on windows
    datapath = 'Y:/research/Experiments/HF-Experiments/GEL-data/';
else
    disp('Platform not supported')
end

end


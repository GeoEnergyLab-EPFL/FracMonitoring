function [datapath] = pathbyarchitecture(servername)
%PATHBYARCHITECTURE builds path depending on computer architecture
%   PATHBYARCHITECTURE returns the HF-experiment server path corresponding
%   to the computer architecture on which the function is running for
%   server SERVERNAME, 
%
% Thomas Blum, Geo-Energy Lab, EPFL, December 2018
%
% NEEDS EXTRA CHECKS FOR WINDOWS PLATFORM SUPPORT

if ismac % path on macos
    [~, username] = system('whoami');
    if strcmp(servername,'enac1files')
        datapath = ['/Users/' username(1:end-1) '/ENACdrives/gel_on_enac1'...
            '/research/Experiments/HF-Experiments/GEL-data/'];
    elseif strcmp(servername,'gel-nas1')
        datapath = ['/Users/' username(1:end-1) '/ENACdrives/'...
            'hf-experiments_on_gel-nas1/'];
    else
        disp('Invalid server')
    end
elseif isunix % path on unix/linus os
    [~, uid] = system('id -u');
    [~, username] = system('whoami');
    rootpath = ['/run/user/' uid(1:end-1) '/gvfs/smb-share:domain=INTRANET,'...
        'server=' servername '.epfl.ch,share='];
    if strcmp(servername,'enac1files')
        datapath = [rootpath 'gel,user=' username(1:end-1) ...
            '/research/Experiments/HF-Experiments/GEL-data/'];
    elseif strcmp(servername,'gel-nas1')
        datapath = [rootpath 'hf-experiments,user=' username(1:end-1) '/'];
    else
        disp('Invalid server')
    end
elseif ispc % path on windows
    datapath = 'X:/';
else
    disp('Platform not supported')
end
end


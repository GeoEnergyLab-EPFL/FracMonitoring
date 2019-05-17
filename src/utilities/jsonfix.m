%% fix for NaN values in JSON files on gel-nas1
% initial parameters
clearvars
close all
home

datapath = pathbyarchitecture('gel-nas1');
FolderInfo = dir(datapath);

% JSON header starts at folder 64
for i_fold = 72:length(FolderInfo)
    tmpInfo = dir([datapath FolderInfo(i_fold).name]);
    for i_file = 1:length(tmpInfo)
        if strfind(tmpInfo(i_file).name,'.json')
            filename = [datapath FolderInfo(i_fold).name '/' tmpInfo(i_file).name];
            % display file being worked on
            disp(filename)
            % extract json header
            jsonhdr = jsondecode(fileread(filename));
            % fix NaN in tubing diameter
            jsonhdr.WellInfos.Tubing.OD = 3.18;
            jsonhdr.WellInfos.Tubing.ID = 1.76;
            % fix passive trigger threshold
            jsonhdr.PassiveAcousticInfos.TriggerThreshold_V_ = 0;
            
            % reencode structure to JSON (NaN is encoded as 'null' by Matlab)
            jsontxt = jsonencode(jsonhdr);
            
            % fix backslashes as special characters before writing to file
            startpos = strfind(jsontxt,'FilePath');
            jsonbegin = jsontxt(1:startpos-1);
            jsonend = jsontxt(startpos:end);
            jsonend = replace(jsonend,'\','\\');
%             backslashpos = regexp(jsonend,'(\w\\\w)');
%             for i_pos = fliplr(backslashpos)
%                 jsonend = [jsonend(1:i_pos) '\' jsonend(i_pos+1:end)];
%             end
            jsontxt2 = [jsonbegin jsonend];
            
            % open new file
            fileID = fopen(filename,'w');
            fprintf(fileID,jsontxt2);
            fclose(fileID);
        end
    end
end
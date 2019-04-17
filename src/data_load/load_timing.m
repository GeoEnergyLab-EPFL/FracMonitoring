function AcqTime = load_timing(filename)
% function to load the low acquisition rate signals from pumps, pressure
% controllers, pressure gauges and other sensors

% check if CSV file exists
fcsv =  strrep(filename,'.bin','.txt');
if ~isfile(fcsv)
        disp('no corresponding txt timing data file found')
    return
end

% read file and save into array
fid = fopen(fcsv);
dataout = textscan(fid,'%{yy-MM-dd}D%{HH:mm:ss}D','Delimiter',',');
fclose(fid);

% create time data
AcqTime = dataout{1};
AcqTime.Format = 'yy-MM-dd HH:mm:ss';
AcqTime = AcqTime + timeofday(dataout{2});

end
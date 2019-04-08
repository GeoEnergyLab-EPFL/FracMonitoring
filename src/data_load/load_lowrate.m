function  [dataout, lowrateTime, varargout] = load_lowrate(filename)
% function to load the low acquisition rate signals from pumps, pressure
% controllers, pressure gauges and other sensors

% check if CSV file exists
fcsv =  strrep(filename,'.bin','.csv');
if ~isfile(fcsv)
        disp('no corresponding low rate CSV data file found')
    return
end

% read first line header
fid = fopen(fcsv);
csvhdr = textscan(fid,'%s',1,'Delimiter',newline);
fclose(fid);

% read data and save into array
dataout = csvread(fcsv,1,0);

% convert header to string 
hdrtmp = strrep(string(csvhdr{:}),',,',', ,');
%varargout(1,:) = strsplit(stmp,',');

% create time vector
tmpTime(:,1) = floor(dataout(:,1)/1E4);
tmpTime(:,2) = floor((dataout(:,1)-tmpTime(:,1)*1E4)/1E2);
tmpTime(:,3) = dataout(:,1)-tmpTime(:,1)*1E4-tmpTime(:,2)*1E2;
lengthTime = size(tmpTime,1);

% look for change in date
I_23 = find(tmpTime(:,1)==23);
I_0 = find(tmpTime(:,1)==0);
i_changetime = I_0(ismember(I_0,I_23+1));

% get start date info from JSON header
jsonhdr = load_header(strrep(filename,'.bin','.json'));
startdate = datetime(jsonhdr.TestInfos.Date_Time,'Format','yyyy-MM-dd;HH:mm:ss');

% start first day 
lowrateTime = datetime(startdate.Year*ones(lengthTime,1),...
    startdate.Month*ones(lengthTime,1),startdate.Day*ones(lengthTime,1),...
    tmpTime(:,1),tmpTime(:,2),tmpTime(:,3));
if i_changetime
    for ii = i_changetime
        lowrateTime(ii:end).Day = lowrateTime(ii:end).Day+1;
    end
end

% optional outputs
nout = max(nargout,1)-2;
varargout = cell(nout,1);
for ii = 1:1
    % only optional output currently is the header line
    varargout(ii) = {strsplit(hdrtmp,',')};
end

end
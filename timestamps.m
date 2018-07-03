function [ AcqTime ] = timestamps(CellTimes,datayear,datamonth,dataday)
%TIMESTAMPS creates a timestamp variable from a timing text file
%   TIMESTAMPS builds the ACQTIME timestamp variable by loading the
%   HH:mm:ss timestamps saves in the ASCII file at FILEPATH and adds the
%   year, month and day from DATAYEAR, DATAMONTH and DATADAY
%
% Thomas Blum, Geo-Energy Lab, EPFL, June 2018

AcqTime = datetime(CellTimes,'Format','HH:mm:ss');  % time in HH:mm:ss
AcqTime.Year = 2000+datayear;
AcqTime.Month = datamonth;
AcqTime.Day = dataday;

end


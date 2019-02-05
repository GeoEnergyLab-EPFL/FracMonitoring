function [transducers,block,json_header]=load_header(filename)
% TODO: write code
% purpose - load json header and create block objects, plattens &
% transducers objects

json_header = jsondecode(fileread(filename));

transducers = [];
block = [250 250 250];

end

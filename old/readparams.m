function s = readparams(filename) %#ok<STOUT>
% READPARAMS   read parameters of scan from file
% READPARAMS(PATH) reads parameter data from the set of measurements in
% PATH (general path without params extension) and returns number of
% samples NS, frequencies NF, traces to stack NT, first frequency FSTART,
% amplification of source AMPL, first angular position STRTPOS, position 
% increment DEGSTEP, and number of positions NBSTEP
%
% Thomas Blum, March 2015

% make sure file exists
if exist([filename '_params.txt'],'file') == 0
    return
end

% variable names to look for
varnames  = {'type','ns','nf','nt','Fs','fstart','ampl','strtpos','step',...
    'nbsteps','degstep','wvlt','dist'};

s = struct;
% read file line by line
fid = fopen([filename '_params.txt'],'r');
tline = fgetl(fid);
while ischar(tline)
    % search through lines and find parameters
    idx =  regexp(tline,'\t');
    if idx == 1
        str = '';
        value  = '';
    elseif idx == length(tline);
        str = tline(1:idx-1);
        value  = '';
    else
        str = tline(1:idx-1);
        value = tline(idx+1:end);
    end
    % fix: change old field "degstep" into more general "step"
    if strcmp('degstep',str)
        str = 'step';
    end
    % write value as string or double depending on type
    if sum(strcmp(str,varnames))
        if prod(isstrprop(value,'alpha'))
            s.(str) = value;
        else
            s.(str) = eval(value);
        end
    end
    % get next line
    tline = fgetl(fid);
end
fclose(fid);

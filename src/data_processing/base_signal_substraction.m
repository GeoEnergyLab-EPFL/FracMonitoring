function [dataSubs] = base_signal_substraction(data,base_sequence)
% subtract a reference signal to the data in order to look at differences
% only
dataSubs = bsxfun(@minus,data,base_sequence);

end

function [dummymode, pp] = initParallelPort(varargin)

% set default parameters
parport = "/dev/parport0";
dummymode = 0;
pp = nan;
% read in the parameters from the function arguments
for i=1:length(varargin)
    switch cell2mat(varargin(i))
        case 'parport'
            parport = varargin{i+1};
        case 'dummymode'
            dummymode = varargin{i+1};
    end
end

if ~dummymode
    % load the instrument control package
    pkg load instrument-control
    
    %% open the port
    pp = parallel(parport, 0);
    pp_data(pp,0)
end

end

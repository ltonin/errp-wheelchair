function [dataout, channel_index] = errp_fill_layout(datain, channel_list, channel_locations)
% ERRP_FILL_LAYOUT  Re-arranges the input data [datain] according to a neurophysiological
% layout given as argument [channel_locations]. 
%
%   [dataout, channel_index] = errp_fill_layout(datain, channel_list, channel_locations [, zerofill])
%   takes the DATAIN [S samples x N channels] and re-arranges the channel positions
%   according to channel names provided with CHANNEL_LIST and present in
%   CHANNEL_LOCATIONS [M channels]. The size of DATAOUT will be [S samples x M channels].
%
%   EXAMPLE:
%       load('chanlocs64.mat'); % Loading 64 channel locations (10-20 system) in 'chalocs'
%       data = randn(100, 32);  % [100 samples x 32 channels]
%       chanlist = {'FP1', 'FP2', 'FZ', 'FC5', 'FC1', 'FC2', 'FC6', 'C3', ...
%                    'CZ', 'C4' , 'CP5', 'CP1', 'CP2', 'CP6', 'P3', 'Pz', ...
%                    'P4', 'F1' , 'F2' , 'FC3', 'FCZ', 'FC4', 'C5', 'C1', ...
%                    'C2', 'C6', 'CP3', 'CP4', 'P5' , 'P1' , 'P2' , 'P6'};
%       data_filled = errp_fill_layout(data, chanlist, chanlocs);
%
%   data_filled is now [100 samples x 64 channels], with the channels
%   re-arranged according to chanlocs64.mat and the non-existent channels
%   set to 0.

    if (length(channel_list) ~= size(datain, 2))
        error('chk:chans', 'Length of channel_list and number of channels in datain must be the same')
    end
    nsamples  = size(datain, 1);
    nchannels = size(datain,2);
    nlocations = length(channel_locations);
    dataout = zeros(nsamples, nlocations);

    chanindex = nan(nchannels, 1);
    for chId = 1:nchannels
        cname = char(channel_list(chId));

        for lId = 1:nlocations
            if(strcmpi(channel_locations(lId).labels, cname) == true)
                chanindex(chId) = lId;
                break;
            end
        end
    end
    
    dataout(:, chanindex) = datain;
    channel_index = chanindex;

end
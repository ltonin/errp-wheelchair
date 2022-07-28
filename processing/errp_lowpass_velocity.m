function out = errp_lowpass_velocity(in, fs, cutoff, order)

    if nargin == 2
        order  = 4;
        cutoff = 1;
    end

    if nargin == 3
        order  = 4;
    end

    %% Interpolate over nan
    out = in;
    if( sum(isnan(in)) > 0)
        s_x = 1:length(in);
        out(isnan(in)) = interp1(s_x(~isnan(in)),in(~isnan(in)),s_x(isnan(in)), 'spline');
    end

    %% Low-pass filter
    [b, a] = butter(order, cutoff*2/fs);
    out = filtfilt(b, a, out);
    
end
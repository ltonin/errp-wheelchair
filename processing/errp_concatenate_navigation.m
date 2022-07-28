function [cmd, pos, vel, events, labels] = errp_concatenate_navigation(files)

    %warning('off', 'backtrace');
    numfiles = length(files);
    
    % Getting size info to allocate memory and speedup the concatenation
    datasize   = get_data_size(files);
    
    nsamples = sum(datasize, 1);

    cmd = nan(nsamples, 3);
    pos = nan(nsamples, 3);
    vel = nan(nsamples, 3);
    Rk  = nan(nsamples, 1);                    % Run
    
    TYP = []; POS = []; DUR = []; ERR = [];
    runId = 1;

    fileseek = 1;
    for fId = 1:numfiles
    
        cfile = files{fId};
        util_disp_progress(fId, numfiles, '        ');
        cdata = load(cfile);

        % Get current position 
        cstart   = fileseek;
        cstop    = cstart + datasize(fId) - 1;
        

        % Concatenate events
        TYP = cat(1, TYP, cdata.events.TYP);
        DUR = cat(1, DUR, cdata.events.DUR);
        ERR = cat(1, ERR, cdata.events.ERR);
        POS = cat(1, POS, cdata.events.POS + fileseek -1);
        

        % Getting cmd, pos, vel
        cmd(cstart:cstop, 1) = cdata.cmd_x;
        cmd(cstart:cstop, 2) = 0;
        cmd(cstart:cstop, 3) = cdata.cmd_z;

        pos(cstart:cstop, 1) = cdata.pos_x;
        pos(cstart:cstop, 2) = cdata.pos_y;
        pos(cstart:cstop, 3) = cdata.pos_z;

        vel(cstart:cstop, 1) = cdata.vel_x;
        vel(cstart:cstop, 2) = 0;
        vel(cstart:cstop, 3) = cdata.vel_z;

        Rk(cstart:cstop) = runId;

        % Update runId
        runId = runId + 1;
        
        % Update the fileseek position
        fileseek = cstop + 1;
        
        
    end
    
    events.TYP = TYP;
    events.POS = POS;
    events.DUR = DUR;
    events.ERR = ERR;

    labels.Rk = Rk;
    
   
   % warning('on', 'backtrace');
end

function dsizes = get_data_size(filepaths)

    nfiles = length(filepaths);
    
    dsizes = zeros(nfiles, 1);
    
    for fId = 1:nfiles
        cfilepath = filepaths{fId};
        cinfo = whos('-file', cfilepath, 't');
        csize = cinfo.size;    
        dsizes(fId) = csize(2);
    end

end




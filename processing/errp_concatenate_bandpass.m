function [F, events, labels, classifiers, settings] = errp_concatenate_bandpass(files)

    %warning('off', 'backtrace');
    numfiles = length(files);
    
    % Getting size info to allocate memory and speedup the concatenation
    datasize   = get_data_size(files);
    NumSamples = sum(datasize(1, :));
    NumChans   = unique(datasize(2, :));
    
    F  = nan(NumSamples, NumChans);
    Rk = nan(NumSamples, 1);                    % Run
    Mk = nan(NumSamples, 1);                    % Modality
    Pk = nan(NumSamples, 1);                    % Protocol
    Dk = nan(NumSamples, 1);                    % Day
    Wk = nan(NumSamples, 1);                    % Week
    Nk = nan(NumSamples, 1);                    % Month
    Sk = nan(NumSamples, 1);
    Dl = [];                                    % Day label
    
    classifiers = [];
    settings = [];
    TYP = []; POS = []; DUR = []; RAC = []; PLY = [];
    runId = 1;
    currday  = 0;
    cweekId  = 0;
    cmonthId = 0;
    lastday  = [];
    lweek    = [];
    lmonth   = [];
    currsubj = '';
    currsubjId = 0;
    
    fileseek = 1;
    for fId = 1:numfiles
    
        cfile = files{fId};
        util_disp_progress(fId, numfiles, '        ');
        cdata = load(cfile);

        % Get current position 
        cstart   = fileseek;
        cstop    = cstart + datasize(1, fId) - 1;
        
        % Get subject
        if(strcmp(cdata.settings.info.subject, currsubj) == false)
            currsubjId = currsubjId + 1;
            currsubj = cdata.settings.info.subject;
        end

        % Get run modality
        if(strcmpi(cdata.settings.modality.name, 'unknown'))
            error(['Unknown modality type: ' cdata.settings.modality.name]);
        else
            [~, cmodality] = intersect(cdata.settings.modality.legend, cdata.settings.modality.name);
        end
        
        % Get run protocol
        if(strcmpi(cdata.settings.protocol.name, 'unknown'))
            error(['Unknown protocol type: ' cdata.settings.protocol.name]);
        else
            [~, cprotocol] = intersect(cdata.settings.protocol.legend, cdata.settings.protocol.name);
        end
        
        % Get day id and label
        if strcmpi(cdata.settings.info.date, lastday) == false
            currday = currday + 1;
            Dl = cat(1, Dl, cdata.settings.info.date);
            lastday = cdata.settings.info.date;
        end

        % Get week id
        cweek = week(datetime(cdata.settings.info.date, 'InputFormat', 'yyyyMMdd'));
        if isequal(cweek, lweek) == false
            cweekId = cweekId +1;
            lweek = cweek;
        end

        % Get month id
        cmonth = month(datetime(cdata.settings.info.date, 'InputFormat', 'yyyyMMdd'));
        if isequal(cmonth, lmonth) == false
            cmonthId = cmonthId +1;
            lmonth = cmonth;
        end

        % Create labels
        Mk(cstart:cstop) = cmodality;
        Pk(cstart:cstop) = cprotocol;
        Rk(cstart:cstop) = runId;
        Dk(cstart:cstop) = currday;
        Wk(cstart:cstop) = cweekId;
        Nk(cstart:cstop) = cmonthId; 
        Sk(cstart:cstop) = currsubjId; 

        % Concatenate events
        TYP = cat(1, TYP, cdata.events.TYP);
        DUR = cat(1, DUR, cdata.events.DUR);
        POS = cat(1, POS, cdata.events.POS + fileseek -1);
        
        tmp = cdata.events;
        if isfield(tmp, 'RAC')
            RAC = cat(1, RAC, cdata.events.RAC);
            PLY = cat(1, PLY, cdata.events.PLY);
        else
            RAC = cat(1, RAC, nan(length(cdata.events.TYP), 1));
            PLY = cat(1, PLY, nan(length(cdata.events.TYP), 1));
        end

        % Getting psd
        F(cstart:cstop, :) = cdata.P;
        
        % Save the current settings, remove number of sample and filename,
        % and compare it with the previous settings. If different, it
        % raises an error.
        csettings = cdata.settings;
        csettings.data.nsamples = nan;
        csettings.data.filename = nan;
        csettings.date          = nan;
        csettings.modality.name = nan;
        csettings.protocol.name = nan;
        csettings.info          = nan;
        csettings.spatial       = nan;
        
        
        if(isempty(settings))
            settings = csettings;
        end
        

        % Getting classifier
        if isfield(cdata, 'classifier')
            classifiers = cat(1, classifiers, cdata.classifier);
        end

        % Update runId
        runId = runId + 1;
        
        % Update the fileseek position
        fileseek = cstop + 1;
        
        
    end
    
    events.TYP = TYP;
    events.POS = POS;
    events.DUR = DUR;
    events.RAC = RAC;
    events.PLY = PLY;
    
    labels.samples.Mk = Mk;
    labels.samples.Pk = Pk;
    labels.samples.Rk = Rk;
    labels.samples.Dk = Dk;
    labels.samples.Wk = Wk;
    labels.samples.Nk = Nk;
    labels.samples.Sk = Sk;
    
    labels.run.Dl = Dl;
   
   % warning('on', 'backtrace');
end

function dsizes = get_data_size(filepaths)

    nfiles = length(filepaths);
    ndimensions = 2;                            % samples x chans
    
    dsizes = zeros(ndimensions, nfiles);
    
    for fId = 1:nfiles
        cfilepath = filepaths{fId};
        cinfo = whos('-file', cfilepath, 'P');
        dsizes(:, fId) = cinfo.size;    
    end

end



